"""Identidade de variante genomica: parsing, normalizacao e validacao de REF.

Uma coordenada sem assembly nao e identidade. `chr17:43093464 A>G` pode ser
GRCh38 ou GRCh37 e os dois apontam para bases diferentes; por isso a assembly e
obrigatoria em toda funcao publica deste modulo.

Formatos aceitos (subset validado, deliberadamente pequeno):
    - VCF-like posicional: `17 43093464 A T` ou `17-43093464-A-T` ou
      `17:43093464:A:T`.
    - SPDI (NCBI Sequence Position Deletion Insertion): `NC_000017.11:43093463:A:T`,
      posicao 0-based conforme a especificacao SPDI.
    - HGVS g. de substituicao simples: `NC_000017.11:g.43093464A>T`.
    - rsID (`rs80357382`) e ClinVar VCV (`VCV000054425`) sao aceitos apenas como
      IDENTIFICADOR opaco para consulta externa, nunca convertidos localmente em
      coordenada.

HGVS c./p., delins complexos, inversoes, duplicacoes e variantes estruturais NAO
sao parseados aqui. HelixScope nao escreve um parser HGVS caseiro: essas notacoes
sao delegadas ao endpoint HGVS do Ensembl VEP, que possui o mapeamento oficial de
transcritos.

Normalizacao: trim de bases comuns a direita e depois a esquerda ate a
representacao minima, seguido de left-alignment quando ha referencia.

A forma minima adotada e a do SPDI/Ensembl, em que um indel puro fica com um
alelo VAZIO (por exemplo `T -> ''` para a delecao de um T). Nao e a forma
ancorada do VCF, que mantem obrigatoriamente uma base em cada alelo
(`CT -> C`). As duas descrevem a mesma alteracao; HelixScope ACEITA a forma
ancorada como entrada e converte-a para a forma minima, porque a forma minima e
unica e por isso serve de identidade e de chave de cache. A forma usada e sempre
declarada no campo `normalization_method`.

Left-shift de repeticoes so ocorre quando ha FASTA de referencia disponivel,
porque deslocar um indel exige ler a sequencia real.

Nota biologica:
    A mesma alteracao de DNA pode ser escrita de varias formas equivalentes
    (por exemplo uma delecao numa regiao repetitiva). Sem normalizacao explicita,
    a mesma variante recebe identidades diferentes e o cache/ClinVar deixam de
    casar. Normalizacao nao altera a biologia; apenas escolhe uma representacao
    canonica.

Nenhuma funcao aqui importa Streamlit.
"""

from __future__ import annotations

import hashlib
import re
from typing import Any, Mapping, Optional

from . import genome_coordinates, genome_fasta, provenance

DNA_LETTERS: frozenset[str] = frozenset("ACGT")

MAX_ALLELE_NT: int = 1_000
"""Tecto de tamanho de alelo aceito neste subset (SNV e indels pequenos)."""

MAX_INPUT_CHARS: int = 512
"""Tecto de comprimento da string de entrada. Evita input abusivo."""

MAX_LEFT_SHIFT_NT: int = 1_000
"""Tecto de deslocamento no left-alignment. Evita varrer um cromossoma inteiro."""

SOURCE_USER: str = "user_input"
SOURCE_VEP: str = "ensembl_vep"
SOURCE_CLINVAR: str = "ncbi_clinvar"

KIND_SNV: str = "SNV"
KIND_INSERTION: str = "insertion"
KIND_DELETION: str = "deletion"
KIND_MNV: str = "MNV"
KIND_INDEL: str = "indel"
KIND_IDENTITY: str = "identity_only"

NORMALIZATION_TRIMMED: str = "TRIMMED_AND_LEFT_ALIGNED"
NORMALIZATION_TRIM_ONLY: str = "TRIMMED_NO_REFERENCE"
NORMALIZATION_NOT_APPLICABLE: str = "NOT_APPLICABLE"

REF_VALIDATED: str = "REF_VALIDATED"
REF_NOT_CHECKED: str = "REF_NOT_CHECKED"
REF_MISMATCH: str = "REFERENCE_MISMATCH"

RSID_PATTERN: re.Pattern[str] = re.compile(r"^rs(\d{1,12})$", re.IGNORECASE)
VCV_PATTERN: re.Pattern[str] = re.compile(r"^VCV(\d{9})(?:\.(\d{1,3}))?$", re.IGNORECASE)
HGVS_G_SUB_PATTERN: re.Pattern[str] = re.compile(
    r"^(?P<seq>[A-Za-z0-9_.]{3,64}):g\.(?P<pos>\d{1,12})(?P<ref>[ACGT]+)>(?P<alt>[ACGT]+)$",
    re.IGNORECASE,
)
SPDI_PATTERN: re.Pattern[str] = re.compile(
    r"^(?P<seq>[A-Za-z0-9_.]{3,64}):(?P<pos>\d{1,12}):(?P<deleted>[ACGT]*|\d{1,6}):(?P<inserted>[ACGT]*)$",
    re.IGNORECASE,
)
CONTIG_PATTERN: re.Pattern[str] = re.compile(r"^[A-Za-z0-9_.\-]{1,64}$")


class VariantInputError(ValueError):
    """Entrada de variante invalida ou nao suportada.

    Attributes:
        category: INVALID_INPUT, UNSUPPORTED, REFERENCE_MISMATCH ou RESOURCE_LIMIT.
    """

    def __init__(self, message: str, category: str = "INVALID_INPUT") -> None:
        super().__init__(message)
        self.category = str(category or "INVALID_INPUT").strip().upper() or "INVALID_INPUT"


def normalize_allele(text: str) -> str:
    """Normaliza um alelo para DNA maiusculo. Vazio representa indel simbolico.

    Args:
        text: Alelo cru. `.`, `-` e `N/A` sao tratados como ausencia de bases.

    Returns:
        String em A/C/G/T, possivelmente vazia.

    Raises:
        VariantInputError: INVALID_INPUT para bases nao DNA; RESOURCE_LIMIT acima
            de MAX_ALLELE_NT.

    Nota biologica:
        U de RNA e convertido para T porque coordenadas genomicas sao de DNA.
        Bases ambiguas (N, R, Y) nao sao aceitas: uma variante precisa de alelos
        definidos para ter identidade.
    """
    seq = str(text or "").strip().upper().replace("U", "T")
    if seq in {".", "-", "N/A", "*"}:
        return ""
    if len(seq) > MAX_ALLELE_NT:
        raise VariantInputError(
            f"Allele exceeds {MAX_ALLELE_NT} nt. Structural variants are not "
            "supported in this subset.",
            "RESOURCE_LIMIT",
        )
    invalid = set(seq) - DNA_LETTERS
    if invalid:
        raise VariantInputError(
            "Alleles must be A/C/G/T (empty allowed for symbolic indel). "
            f"Refused characters: {''.join(sorted(invalid))}."
        )
    return seq


def classify_variant_kind(ref: str, alt: str) -> str:
    """Classifica a variante por comprimento de alelo. Nao prediz efeito.

    Args:
        ref: Alelo de referencia normalizado.
        alt: Alelo alternativo normalizado.

    Returns:
        SNV, MNV, insertion, deletion ou indel.

    Raises:
        Nenhum.

    Nota biologica:
        A classe e puramente descritiva do comprimento. `deletion` nao implica
        perda de funcao e `SNV` nao implica troca de aminoacido; consequencia
        proteica exige anotacao de transcrito.
    """
    r = str(ref or "")
    a = str(alt or "")
    if len(r) == 1 and len(a) == 1:
        return KIND_SNV
    if len(r) == len(a) and len(r) > 1:
        return KIND_MNV
    if not r and a:
        return KIND_INSERTION
    if r and not a:
        return KIND_DELETION
    if len(a) > len(r):
        return KIND_INSERTION
    if len(r) > len(a):
        return KIND_DELETION
    return KIND_INDEL


def trim_alleles(*, position_0based: int, ref: str, alt: str) -> dict:
    """Reduz ref/alt a representacao minima (forma SPDI, alelo vazio permitido).

    Args:
        position_0based: Posicao 0-based da primeira base de `ref`.
        ref: Alelo de referencia normalizado.
        alt: Alelo alternativo normalizado.

    Returns:
        Dict position_0based, ref, alt, trimmed_right, trimmed_left.

    Raises:
        VariantInputError: INVALID_INPUT se ref e alt forem iguais ou a posicao
            for negativa.

    Nota biologica:
        `CTT -> CT`, `TT -> T` e `T -> ''` descrevem a mesma delecao de um unico
        T. Esta funcao devolve a ultima forma, a minima, porque e a unica
        representacao unica e por isso a unica que serve de identidade. Manter
        bases que nao mudam faria a mesma variante receber identidades
        diferentes conforme quem a escreveu.
    """
    try:
        pos = int(position_0based)
    except (TypeError, ValueError) as exc:
        raise VariantInputError("position_0based must be an integer.") from exc
    if pos < 0:
        raise VariantInputError("position_0based must be >= 0.")
    r = str(ref or "")
    a = str(alt or "")
    if r == a:
        raise VariantInputError(
            "Ref and alt are identical; this is not a variant."
        )
    trimmed_right = 0
    while r and a and len(r) > 0 and len(a) > 0 and r[-1] == a[-1] and (len(r) > 1 or len(a) > 1):
        r = r[:-1]
        a = a[:-1]
        trimmed_right += 1
    trimmed_left = 0
    while r and a and r[0] == a[0] and (len(r) > 1 or len(a) > 1):
        r = r[1:]
        a = a[1:]
        pos += 1
        trimmed_left += 1
    return {
        "position_0based": pos,
        "ref": r,
        "alt": a,
        "trimmed_right": trimmed_right,
        "trimmed_left": trimmed_left,
    }


def _fai_row(fasta_path: str, fai_path: str, contig: str) -> Mapping[str, Any]:
    index = genome_fasta.load_fai(fai_path)
    row = (index.get("contigs") or {}).get(str(contig))
    if not row:
        available = sorted((index.get("contigs") or {}).keys())[:8]
        raise VariantInputError(
            f"Contig {contig} is not present in the reference index. "
            f"Indexed contigs include: {', '.join(available) or 'none'}. "
            "HelixScope does not silently map chr17 to NC_000017.11.",
            "INVALID_INPUT",
        )
    return row


def left_align_indel(
    *,
    fasta_path: str,
    fai_path: str,
    contig: str,
    position_0based: int,
    ref: str,
    alt: str,
) -> dict:
    """Desloca um indel para a posicao mais a esquerda equivalente.

    Args:
        fasta_path: FASTA indexado da assembly.
        fai_path: Indice FAI correspondente.
        contig: Identificador do contig exatamente como no FASTA.
        position_0based: Posicao 0-based apos trim.
        ref: Alelo de referencia apos trim.
        alt: Alelo alternativo apos trim.

    Returns:
        Dict position_0based, ref, alt, shifted_nt.

    Raises:
        VariantInputError: INVALID_INPUT se o contig nao existir; RESOURCE_LIMIT
            se o deslocamento exceder MAX_LEFT_SHIFT_NT.

    Nota biologica:
        Numa regiao homopolimerica (por exemplo AAAAAA) a delecao de um A pode ser
        escrita em varias posicoes. A convencao VCF e a posicao mais a esquerda.
        Sem o FASTA nao e possivel deslocar honestamente, por isso esta funcao
        exige a referencia.
    """
    r = str(ref or "")
    a = str(alt or "")
    if r and a:
        return {"position_0based": int(position_0based), "ref": r, "alt": a, "shifted_nt": 0}
    row = _fai_row(fasta_path, fai_path, contig)
    pos = int(position_0based)
    shifted = 0
    while pos > 0 and shifted < MAX_LEFT_SHIFT_NT:
        previous = genome_fasta.fetch_sequence(fasta_path, row, pos - 1, pos)
        if not previous:
            break
        base = previous.upper()
        candidate_ref = base + r
        candidate_alt = base + a
        if candidate_ref[-1] != candidate_alt[-1]:
            break
        trimmed = trim_alleles(
            position_0based=pos - 1, ref=candidate_ref, alt=candidate_alt
        )
        if trimmed["position_0based"] >= pos:
            break
        pos = int(trimmed["position_0based"])
        r = str(trimmed["ref"])
        a = str(trimmed["alt"])
        shifted += 1
    if shifted >= MAX_LEFT_SHIFT_NT:
        raise VariantInputError(
            f"Left-alignment exceeded {MAX_LEFT_SHIFT_NT} nt. Refusing to scan "
            "further; the input may be a structural variant.",
            "RESOURCE_LIMIT",
        )
    return {"position_0based": pos, "ref": r, "alt": a, "shifted_nt": shifted}


def validate_ref_against_reference(
    *,
    fasta_path: str,
    fai_path: str,
    contig: str,
    position_0based: int,
    ref: str,
) -> dict:
    """Confirma que o alelo REF casa com o FASTA da assembly.

    Args:
        fasta_path: FASTA indexado.
        fai_path: Indice FAI.
        contig: Contig exatamente como no FASTA.
        position_0based: Posicao 0-based da primeira base de REF.
        ref: Alelo de referencia (vazio para insercao pura).

    Returns:
        Dict status, observed, expected, contig_length.

    Raises:
        VariantInputError: INVALID_INPUT se o contig nao existir ou a posicao cair
            fora do contig.

    Nota biologica:
        REF errado e o erro silencioso mais comum em analise de variantes: em
        geral significa que a coordenada e de outra assembly (GRCh37 vs GRCh38).
        Por isso a divergencia devolve REFERENCE_MISMATCH em vez de continuar.
    """
    row = _fai_row(fasta_path, fai_path, contig)
    length = int(row["length"])
    pos = int(position_0based)
    expected = str(ref or "").upper()
    if pos < 0 or pos > length:
        raise VariantInputError(
            f"Position {pos} (0-based) is outside contig {contig} of length {length}.",
            "INVALID_INPUT",
        )
    if not expected:
        return {
            "status": REF_VALIDATED,
            "observed": "",
            "expected": "",
            "contig_length": length,
            "reason": (
                "Pure insertion has no REF bases to compare; the anchor position "
                "was checked against the contig bounds only."
            ),
        }
    if pos + len(expected) > length:
        raise VariantInputError(
            f"REF of {len(expected)} nt at 0-based {pos} runs past the end of "
            f"contig {contig} (length {length}).",
            "INVALID_INPUT",
        )
    observed = genome_fasta.fetch_sequence(
        fasta_path, row, pos, pos + len(expected)
    ).upper()
    if observed != expected:
        return {
            "status": REF_MISMATCH,
            "observed": observed,
            "expected": expected,
            "contig_length": length,
            "reason": (
                f"Reference FASTA has {observed} at {contig} 1-based "
                f"{pos + 1}, but the input declares REF {expected}. This is "
                "usually a wrong assembly (for example GRCh37 coordinates "
                "submitted as GRCh38). HelixScope refuses the variant instead "
                "of annotating the wrong base."
            ),
        }
    return {
        "status": REF_VALIDATED,
        "observed": observed,
        "expected": expected,
        "contig_length": length,
        "reason": "REF matches the indexed reference FASTA at this coordinate.",
    }


def parse_variant_input(text: str) -> dict:
    """Reconhece o subset de formatos suportado. Nao consulta a rede.

    Args:
        text: String de entrada do utilizador.

    Returns:
        Dict com `format` e os campos reconhecidos. Para rsID/VCV devolve
        `format` identifier e `identifier`, sem coordenada.

    Raises:
        VariantInputError: INVALID_INPUT para vazio/comprimento; UNSUPPORTED para
            notacoes que exigem um motor de transcritos (HGVS c./p., delins).

    Nota biologica:
        rsID e VCV nao carregam coordenada por si: o mesmo rsID pode ter
        posicoes diferentes em assemblies diferentes. Por isso sao devolvidos
        como identificador opaco para resolucao por servico oficial.
    """
    raw = str(text or "").strip()
    if not raw:
        raise VariantInputError("Variant input is empty.")
    if len(raw) > MAX_INPUT_CHARS:
        raise VariantInputError(
            f"Variant input exceeds {MAX_INPUT_CHARS} characters.", "RESOURCE_LIMIT"
        )
    rs = RSID_PATTERN.match(raw)
    if rs:
        return {
            "format": "rsid",
            "identifier": f"rs{int(rs.group(1))}",
            "needs_external_resolution": True,
            "reason": (
                "An rsID has no assembly-specific coordinate on its own. It must "
                "be resolved by an official service before annotation."
            ),
        }
    vcv = VCV_PATTERN.match(raw)
    if vcv:
        version = vcv.group(2)
        return {
            "format": "clinvar_vcv",
            "identifier": f"VCV{vcv.group(1)}" + (f".{version}" if version else ""),
            "needs_external_resolution": True,
            "reason": (
                "A ClinVar VCV accession is an identifier, not a coordinate. "
                "HelixScope retrieves the record instead of guessing a position."
            ),
        }
    hgvs = HGVS_G_SUB_PATTERN.match(raw)
    if hgvs:
        ref = normalize_allele(hgvs.group("ref"))
        alt = normalize_allele(hgvs.group("alt"))
        if len(ref) != len(alt):
            raise VariantInputError(
                "Only simple HGVS g. substitutions are parsed locally. Use the "
                "Ensembl VEP HGVS endpoint for delins, dup, inv and del.",
                "UNSUPPORTED",
            )
        return {
            "format": "hgvs_g_substitution",
            "sequence_id": hgvs.group("seq"),
            "position_0based": int(hgvs.group("pos")) - 1,
            "ref": ref,
            "alt": alt,
            "needs_external_resolution": False,
            "reason": "HGVS g. substitution: 1-based position converted to internal 0-based.",
        }
    lowered = raw.lower()
    if ":c." in lowered or ":p." in lowered or ":n." in lowered or ":m." in lowered:
        raise VariantInputError(
            "HGVS c./p./n./m. notation needs an official transcript model. "
            "HelixScope does not implement a house HGVS parser; submit this "
            "notation to the Ensembl VEP HGVS endpoint instead.",
            "UNSUPPORTED",
        )
    spdi = SPDI_PATTERN.match(raw)
    if spdi and ":" in raw and raw.count(":") == 3:
        deleted = spdi.group("deleted")
        if deleted.isdigit():
            raise VariantInputError(
                "SPDI with a numeric deletion length requires the reference "
                "sequence to expand the deleted bases. Provide explicit REF "
                "bases or load the matching assembly.",
                "UNSUPPORTED",
            )
        return {
            "format": "spdi",
            "sequence_id": spdi.group("seq"),
            "position_0based": int(spdi.group("pos")),
            "ref": normalize_allele(deleted),
            "alt": normalize_allele(spdi.group("inserted")),
            "needs_external_resolution": False,
            "reason": "SPDI position is already 0-based per the NCBI specification.",
        }
    # Split on whitespace/colon/comma/pipe first so a standalone '-' or '.'
    # stays readable as an empty allele. Only fall back to '-' as a separator
    # for the compact `17-43093557-C-G` form.
    fields = [item for item in re.split(r"[\s:,|]+", raw) if item]
    if len(fields) != 4:
        fields = [item for item in re.split(r"[\s:\-,|]+", raw) if item]
    if len(fields) == 4:
        contig, position, ref, alt = fields
        if not CONTIG_PATTERN.match(contig):
            raise VariantInputError(
                "Contig token contains unsupported characters."
            )
        if not position.isdigit():
            raise VariantInputError("Position must be a 1-based integer in VCF-like input.")
        return {
            "format": "vcf_like",
            "sequence_id": contig,
            "position_0based": int(position) - 1,
            "ref": normalize_allele(ref),
            "alt": normalize_allele(alt),
            "needs_external_resolution": False,
            "reason": (
                "VCF-like input is 1-based inclusive; converted to internal "
                "0-based half-open."
            ),
        }
    raise VariantInputError(
        "Unrecognized variant format. Supported: VCF-like "
        "'17 43093464 A T', SPDI 'NC_000017.11:43093463:A:T', HGVS g. "
        "substitution 'NC_000017.11:g.43093464A>T', rsID, or ClinVar VCV.",
        "UNSUPPORTED",
    )


def variant_identity_hash(
    *,
    assembly: str,
    contig: str,
    position_0based: int,
    ref: str,
    alt: str,
) -> str:
    """Hash estavel da identidade normalizada, incluindo a assembly.

    Args:
        assembly: Nome/accession da assembly.
        contig: Contig.
        position_0based: Posicao 0-based normalizada.
        ref: REF normalizado.
        alt: ALT normalizado.

    Returns:
        SHA-256 hexadecimal.

    Raises:
        Nenhum.

    Nota biologica:
        A assembly entra no hash de proposito: a mesma coordenada em GRCh37 e
        GRCh38 e uma variante diferente e nunca deve reutilizar cache.
    """
    payload = "|".join(
        [
            str(assembly or "").strip().upper(),
            str(contig or "").strip(),
            str(int(position_0based)),
            str(ref or "").upper(),
            str(alt or "").upper(),
        ]
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def build_variant(
    *,
    text: str,
    assembly: str,
    accession: str = "",
    fasta_path: str = "",
    fai_path: str = "",
    assembly_sha256: str = "",
    source: str = SOURCE_USER,
    contig_override: str = "",
) -> dict:
    """Constroi a identidade completa de uma variante. Assembly e obrigatoria.

    Args:
        text: Entrada crua num dos formatos suportados.
        assembly: Nome da assembly declarada (ex. GRCh38.p14). Obrigatorio.
        accession: Accession da assembly quando conhecido (ex. GCF_000001405.40).
        fasta_path: FASTA READY para validar REF e alinhar indels (opcional).
        fai_path: Indice FAI correspondente (opcional).
        assembly_sha256: SHA-256 do FASTA READY quando existir.
        source: Origem do registo (user_input, ensembl_vep, ncbi_clinvar).
        contig_override: Contig a usar quando o input traz um alias que o
            utilizador confirmou explicitamente.

    Returns:
        Dict imutavel em conteudo com identidade, normalizacao, validacao de REF
        e proveniencia. `effect` permanece None: este modulo nao prediz efeito.

    Raises:
        VariantInputError: INVALID_INPUT sem assembly, formato invalido, ou
            REFERENCE_MISMATCH quando o FASTA contradiz o REF declarado.

    Nota biologica:
        O objeto devolvido descreve apenas ONDE e QUAL e a troca de bases.
        Consequencia molecular, patogenicidade e efeito proteico vem de motores
        externos (Ensembl VEP) e de evidencia submetida (ClinVar), nunca daqui.
    """
    declared = str(assembly or "").strip()
    if not declared:
        raise VariantInputError(
            "Assembly is required. A genomic coordinate without an assembly is "
            "not an identity: chr17:43093464 A>G means different bases in "
            "GRCh37 and GRCh38.",
            "INVALID_INPUT",
        )
    parsed = parse_variant_input(text)
    if parsed.get("needs_external_resolution"):
        return {
            "kind": KIND_IDENTITY,
            "status": "IDENTIFIER_ONLY",
            "format": parsed["format"],
            "identifier": parsed["identifier"],
            "assembly": declared,
            "accession": str(accession or ""),
            "contig": "",
            "position_0based": None,
            "position_1based": None,
            "ref": "",
            "alt": "",
            "variant_kind": None,
            "normalization": NORMALIZATION_NOT_APPLICABLE,
            "ref_validation": REF_NOT_CHECKED,
            "ref_validation_reason": str(parsed.get("reason") or ""),
            "identity_hash": "",
            "effect": None,
            "effect_status": "NOT_COMPUTED",
            "effect_note": (
                "HelixScope does not predict variant effect. Consequence terms "
                "come from Ensembl VEP; clinical assertions come from ClinVar "
                "submitters."
            ),
            "source": str(source or SOURCE_USER),
            "assembly_sha256": str(assembly_sha256 or ""),
            "internal_convention": genome_coordinates.INTERNAL_CONVENTION,
            "display_convention": genome_coordinates.DISPLAY_CONVENTION,
            "retrieved_at_utc": provenance.utc_now(),
            "software_version": provenance.HELIXSCOPE_VERSION,
        }
    contig = str(contig_override or parsed.get("sequence_id") or "").strip()
    if not CONTIG_PATTERN.match(contig):
        raise VariantInputError("Contig identifier is missing or not a safe token.")
    ref = str(parsed.get("ref") or "")
    alt = str(parsed.get("alt") or "")
    if not ref and not alt:
        raise VariantInputError("Ref and alt cannot both be empty.")
    trimmed = trim_alleles(
        position_0based=int(parsed["position_0based"]), ref=ref, alt=alt
    )
    position = int(trimmed["position_0based"])
    ref = str(trimmed["ref"])
    alt = str(trimmed["alt"])
    normalization = NORMALIZATION_TRIM_ONLY
    shifted = 0
    ref_status = REF_NOT_CHECKED
    ref_reason = (
        "No READY reference FASTA was supplied, so REF was not verified and "
        "indels were not left-aligned. Identity is trim-normalized only."
    )
    ref_observed = ""
    has_reference = bool(fasta_path and fai_path)
    if has_reference:
        aligned = left_align_indel(
            fasta_path=fasta_path,
            fai_path=fai_path,
            contig=contig,
            position_0based=position,
            ref=ref,
            alt=alt,
        )
        position = int(aligned["position_0based"])
        ref = str(aligned["ref"])
        alt = str(aligned["alt"])
        shifted = int(aligned["shifted_nt"])
        normalization = NORMALIZATION_TRIMMED
        validation = validate_ref_against_reference(
            fasta_path=fasta_path,
            fai_path=fai_path,
            contig=contig,
            position_0based=position,
            ref=ref,
        )
        ref_status = str(validation["status"])
        ref_reason = str(validation["reason"])
        ref_observed = str(validation["observed"])
        if ref_status == REF_MISMATCH:
            raise VariantInputError(ref_reason, "REFERENCE_MISMATCH")
    display = genome_coordinates.internal_to_display(
        position, position + max(len(ref), 1)
    )
    return {
        "kind": KIND_IDENTITY,
        "status": "NORMALIZED",
        "format": parsed["format"],
        "identifier": "",
        "assembly": declared,
        "accession": str(accession or ""),
        "contig": contig,
        "position_0based": position,
        "position_1based": display.get("start_1based"),
        "ref": ref,
        "alt": alt,
        "variant_kind": classify_variant_kind(ref, alt),
        "normalization": normalization,
        "normalization_method": (
            "Right then left trim of common bases down to the minimal SPDI-style "
            "representation (a pure indel keeps one empty allele), followed by "
            "left-alignment against the reference FASTA when one is available. "
            "VCF anchored input is accepted and converted to this form."
        ),
        "trimmed_left": int(trimmed["trimmed_left"]),
        "trimmed_right": int(trimmed["trimmed_right"]),
        "left_shifted_nt": shifted,
        "ref_validation": ref_status,
        "ref_validation_reason": ref_reason,
        "ref_observed": ref_observed,
        "identity_hash": variant_identity_hash(
            assembly=declared,
            contig=contig,
            position_0based=position,
            ref=ref,
            alt=alt,
        ),
        "effect": None,
        "effect_status": "NOT_COMPUTED",
        "effect_note": (
            "HelixScope does not predict variant effect. Consequence terms come "
            "from Ensembl VEP; clinical assertions come from ClinVar submitters."
        ),
        "source": str(source or SOURCE_USER),
        "assembly_sha256": str(assembly_sha256 or ""),
        "internal_convention": genome_coordinates.INTERNAL_CONVENTION,
        "display_convention": genome_coordinates.DISPLAY_CONVENTION,
        "retrieved_at_utc": provenance.utc_now(),
        "software_version": provenance.HELIXSCOPE_VERSION,
    }


def variant_region_string(variant: Mapping[str, Any]) -> str:
    """Formata a variante como regiao Ensembl VEP `chr:start-end:strand`.

    Args:
        variant: Objeto de `build_variant` com status NORMALIZED.

    Returns:
        String de regiao 1-based aceita pelo endpoint region do VEP.

    Raises:
        VariantInputError: INVALID_INPUT se a variante nao estiver normalizada.

    Nota biologica:
        O VEP usa coordenadas 1-based inclusive. Para insercao pura, a spec do
        endpoint pede start = end + 1, indicando o ponto de insercao entre bases.
    """
    if str(variant.get("status") or "") != "NORMALIZED":
        raise VariantInputError(
            "Only a normalized variant with an explicit coordinate can build a "
            "VEP region string."
        )
    contig = str(variant.get("contig") or "")
    start_1 = int(variant.get("position_1based") or 0)
    ref = str(variant.get("ref") or "")
    if not ref:
        return f"{contig}:{start_1}-{start_1 - 1}:1"
    end_1 = start_1 + len(ref) - 1
    return f"{contig}:{start_1}-{end_1}:1"


def ensembl_contig_name(contig: str) -> str:
    """Devolve o nome de contig no estilo Ensembl, sem inventar equivalencias.

    Args:
        contig: Contig como escrito pelo utilizador ou pelo FASTA.

    Returns:
        Nome sem prefixo `chr` quando o resto e um cromossoma humano canonico
        (1-22, X, Y, MT/M). Caso contrario devolve o texto original.

    Raises:
        Nenhum.

    Nota biologica:
        Ensembl nomeia cromossomas humanos como `17`; UCSC usa `chr17`; RefSeq
        usa `NC_000017.11`. Remover `chr` de um cromossoma canonico e uma
        equivalencia documentada e segura. Um accession RefSeq NAO e convertido
        aqui: essa traducao exige a tabela oficial da assembly.
    """
    text = str(contig or "").strip()
    if not text:
        return ""
    bare = text[3:] if text.lower().startswith("chr") else text
    canonical = {str(index) for index in range(1, 23)} | {"X", "Y", "MT", "M"}
    if bare.upper() in canonical:
        return "MT" if bare.upper() in {"MT", "M"} else bare.upper() if bare.upper() in {"X", "Y"} else bare
    return text
