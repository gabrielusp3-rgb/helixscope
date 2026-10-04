"""Diagnostico OpenCL para Cas-OFFinder. Nao instala runtimes.

Cas-OFFinder 2.x (Park e Bae, 2014) usa OpenCL. No Windows o modo C (CPU)
ainda exige um ICD OpenCL funcional. Este modulo so observa: loader, registo
Khronos, ICDs oficiais conhecidos, plataformas via clGetPlatformIDs.

Nao descarrega SDKs, nao copia OpenCL.dll de foruns, nao altera o System32.
Nenhuma funcao importa Streamlit.

Nota biologica:
    A ausencia de dispositivo OpenCL impede a busca Cas-OFFinder. Nao e um
    resultado de off-targets vazio e nao autoriza rotulo genome-wide.
"""

from __future__ import annotations

import ctypes
import os
import platform
import subprocess
from typing import Any, Mapping, Optional

from . import provenance

OFFICIAL_INTEL_CPU_RUNTIME_NOTE: str = (
    "Intel CPU Runtime for OpenCL Applications is the vendor CPU ICD for "
    "Intel Core/Xeon. Official documentation (retrieved 2026-08-28): "
    "https://www.intel.com/content/www/us/en/developer/articles/technical/"
    "intel-cpu-runtime-for-opencl-applications-with-sycl-support.html "
    "Windows installer example: w_opencl_runtime_p_2026.0.0.946.exe. "
    "Requires administrator install. HelixScope does not run that installer."
)

OFFICIAL_INTEL_GPU_NOTE: str = (
    "Intel Graphics OpenCL ICD ships with the Intel Graphics Driver. From "
    "the February 2020 graphics packages onward, CPU OpenCL is no longer "
    "bundled in that driver (Intel OpenCL Runtimes for Intel Processors). "
    "An unregistered GPU ICD in DriverStore is not a working CPU device."
)

# ICD de graficos Intel observado nesta maquina: apontar OCL_ICD_FILENAMES
# para este ficheiro fez Cas-OFFinder 2.4.1 bloquear sem stdout (C e G).
HANGING_ICD_BASENAMES: frozenset[str] = frozenset(
    {
        "intel_opencl_icd64.dll",
        "igfx11cmrt64.dll",
    }
)

CPU_ICD_BASENAMES: frozenset[str] = frozenset(
    {
        "intelocl64.dll",
        "intelocl32.dll",
    }
)


class OpenCLError(ValueError):
    """Falha de diagnostico OpenCL.

    Attributes:
        category: INVALID_INPUT ou ERROR.
    """

    def __init__(self, message: str, category: str = "ERROR") -> None:
        super().__init__(message)
        self.category = str(category or "ERROR").strip().upper() or "ERROR"


def hardware_summary() -> dict:
    """CPU/GPU minimos para escolher o runtime oficial.

    Args:
        Nenhum.

    Returns:
        Dict architecture, processor, gpu_names, os.

    Raises:
        Nenhum.
    """
    gpus: list[str] = []
    cpu_name = platform.processor() or ""
    cpu_cores = os.cpu_count()
    if os.name == "nt" and os.environ.get("HELIXSCOPE_OPENCL_PROBE_WMI") == "1":
        gpus = _windows_gpu_names()
        wmi_cpu = _windows_cpu_name()
        if wmi_cpu:
            cpu_name = wmi_cpu
    return {
        "os": platform.system(),
        "os_release": platform.release(),
        "architecture": platform.machine(),
        "processor": cpu_name,
        "logical_cpus": cpu_cores,
        "gpu_names": gpus,
        "note": (
            "GPU is not required. Cas-OFFinder device C is preferred when a "
            "CPU OpenCL ICD is registered. Absence of a discrete GPU is not "
            "an artificial blocker."
        ),
    }


def _windows_cpu_name() -> str:
    try:
        completed = subprocess.run(
            [
                "powershell",
                "-NoProfile",
                "-Command",
                "(Get-CimInstance Win32_Processor).Name",
            ],
            capture_output=True,
            text=True,
            timeout=8,
            check=False,
            shell=False,
        )
    except (OSError, subprocess.TimeoutExpired):
        return ""
    return (completed.stdout or "").strip().splitlines()[0] if completed.stdout else ""


def _windows_gpu_names() -> list[str]:
    try:
        completed = subprocess.run(
            [
                "powershell",
                "-NoProfile",
                "-Command",
                "Get-CimInstance Win32_VideoController | ForEach-Object { $_.Name }",
            ],
            capture_output=True,
            text=True,
            timeout=8,
            check=False,
            shell=False,
        )
    except (OSError, subprocess.TimeoutExpired):
        return []
    names = []
    for line in (completed.stdout or "").splitlines():
        text = line.strip()
        if text:
            names.append(text[:120])
    return names


def icd_loader_path() -> str:
    """Caminho do loader OpenCL do sistema, se existir.

    Args:
        Nenhum.

    Returns:
        Caminho ou string vazia.

    Raises:
        Nenhum.
    """
    if os.name == "nt":
        system = os.environ.get("SystemRoot") or r"C:\Windows"
        candidate = os.path.join(system, "System32", "OpenCL.dll")
        if os.path.isfile(candidate):
            return candidate
        lower = os.path.join(system, "System32", "opencl.dll")
        if os.path.isfile(lower):
            return lower
        return ""
    for path in (
        "/usr/lib/x86_64-linux-gnu/libOpenCL.so.1",
        "/usr/lib/libOpenCL.so.1",
        "/usr/local/lib/libOpenCL.so.1",
    ):
        if os.path.isfile(path):
            return path
    return ""


def khronos_vendor_icds() -> list[dict]:
    """ICDs registados no mecanismo oficial da plataforma.

    Args:
        Nenhum.

    Returns:
        Lista de dicts path, registered, exists.

    Raises:
        Nenhum.
    """
    rows: list[dict] = []
    if os.name == "nt":
        rows.extend(_windows_khronos_vendors())
    else:
        vendor_dir = "/etc/OpenCL/vendors"
        if os.path.isdir(vendor_dir):
            try:
                names = os.listdir(vendor_dir)
            except OSError:
                names = []
            for name in names:
                if not name.endswith(".icd"):
                    continue
                icd_file = os.path.join(vendor_dir, name)
                dll = ""
                try:
                    with open(icd_file, "r", encoding="utf-8", errors="replace") as handle:
                        dll = handle.readline().strip()
                except OSError:
                    dll = ""
                rows.append(
                    {
                        "registry_entry": name,
                        "path": dll,
                        "registered": True,
                        "exists": bool(dll) and os.path.isfile(dll),
                    }
                )
    return rows


def _windows_khronos_vendors() -> list[dict]:
    try:
        import winreg
    except ImportError:
        return []
    rows: list[dict] = []
    for hive, hive_name in (
        (getattr(winreg, "HKEY_LOCAL_MACHINE", None), "HKLM"),
        (getattr(winreg, "HKEY_CURRENT_USER", None), "HKCU"),
    ):
        if hive is None:
            continue
        for key_path in (
            r"SOFTWARE\Khronos\OpenCL\Vendors",
            r"SOFTWARE\WOW6432Node\Khronos\OpenCL\Vendors",
        ):
            try:
                key = winreg.OpenKey(hive, key_path)
            except OSError:
                continue
            try:
                index = 0
                while True:
                    try:
                        name, _value, _typ = winreg.EnumValue(key, index)
                    except OSError:
                        break
                    rows.append(
                        {
                            "registry_entry": f"{hive_name}\\{key_path}\\{name}",
                            "path": str(name),
                            "registered": True,
                            "exists": os.path.isfile(str(name)),
                        }
                    )
                    index += 1
            finally:
                winreg.CloseKey(key)
    return rows


def enumerate_platforms() -> dict:
    """Chama clGetPlatformIDs no loader do sistema (ICDs registados).

    Args:
        Nenhum.

    Returns:
        Dict n_platforms, loader_loadable, devices (lista curta).

    Raises:
        Nenhum.
    """
    loader = icd_loader_path()
    if not loader:
        return {
            "loader_loadable": False,
            "n_platforms": 0,
            "platforms": [],
            "reason": "OpenCL ICD loader library was not found.",
        }
    try:
        if os.name == "nt":
            ocl = ctypes.WinDLL(loader)
        else:
            ocl = ctypes.CDLL(loader)
    except OSError as exc:
        return {
            "loader_loadable": False,
            "n_platforms": 0,
            "platforms": [],
            "reason": f"OpenCL loader could not be loaded: {exc}.",
        }
    clGetPlatformIDs = ocl.clGetPlatformIDs
    clGetPlatformIDs.restype = ctypes.c_int32
    clGetPlatformIDs.argtypes = [
        ctypes.c_uint32,
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_uint32),
    ]
    count = ctypes.c_uint32(0)
    err = int(clGetPlatformIDs(0, None, ctypes.byref(count)))
    n_plat = int(count.value)
    platforms: list[dict] = []
    if err == 0 and n_plat > 0:
        plat_array = (ctypes.c_void_p * n_plat)()
        err2 = int(clGetPlatformIDs(n_plat, plat_array, None))
        if err2 == 0:
            for handle in plat_array:
                platforms.append(_platform_info(ocl, int(handle or 0)))
    return {
        "loader_loadable": True,
        "loader_path_basename": os.path.basename(loader),
        "n_platforms": n_plat,
        "clGetPlatformIDs_error": err,
        "platforms": platforms,
        "reason": (
            ""
            if n_plat
            else "ICD loader is present but reports zero OpenCL platforms. "
            "No vendor ICD is registered."
        ),
    }


def _platform_info(ocl: Any, handle: int) -> dict:
    if not handle:
        return {"vendor": "", "name": "", "version": "", "devices": []}
    name = _cl_str(ocl, handle, 0x0902)  # CL_PLATFORM_NAME
    vendor = _cl_str(ocl, handle, 0x0903)  # CL_PLATFORM_VENDOR
    version = _cl_str(ocl, handle, 0x0901)  # CL_PLATFORM_VERSION
    devices = _devices_for_platform(ocl, handle)
    return {
        "name": name[:80],
        "vendor": vendor[:80],
        "version": version[:80],
        "devices": devices,
    }


def _cl_str(ocl: Any, handle: int, param: int) -> str:
    fn = ocl.clGetPlatformInfo
    fn.restype = ctypes.c_int32
    fn.argtypes = [
        ctypes.c_void_p,
        ctypes.c_uint32,
        ctypes.c_size_t,
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_size_t),
    ]
    size = ctypes.c_size_t(0)
    err = int(fn(handle, param, 0, None, ctypes.byref(size)))
    if err != 0 or int(size.value) <= 1:
        return ""
    buf = ctypes.create_string_buffer(int(size.value))
    err = int(fn(handle, param, int(size.value), buf, None))
    if err != 0:
        return ""
    return buf.value.decode("utf-8", errors="replace").strip()


def _devices_for_platform(ocl: Any, platform_handle: int) -> list[dict]:
    clGetDeviceIDs = ocl.clGetDeviceIDs
    clGetDeviceIDs.restype = ctypes.c_int32
    clGetDeviceIDs.argtypes = [
        ctypes.c_void_p,
        ctypes.c_uint64,
        ctypes.c_uint32,
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_uint32),
    ]
    n_dev = ctypes.c_uint32(0)
    # CL_DEVICE_TYPE_ALL = 0xFFFFFFFF
    err = int(
        clGetDeviceIDs(platform_handle, 0xFFFFFFFF, 0, None, ctypes.byref(n_dev))
    )
    if err != 0 or int(n_dev.value) == 0:
        return []
    count = int(n_dev.value)
    dev_array = (ctypes.c_void_p * count)()
    err = int(clGetDeviceIDs(platform_handle, 0xFFFFFFFF, count, dev_array, None))
    if err != 0:
        return []
    rows: list[dict] = []
    for handle in dev_array:
        rows.append(_device_info(ocl, int(handle or 0)))
    return rows


def _device_info(ocl: Any, handle: int) -> dict:
    if not handle:
        return {"name": "", "type": "", "vendor": "", "version": ""}
    name = _cl_device_str(ocl, handle, 0x102B)  # CL_DEVICE_NAME
    vendor = _cl_device_str(ocl, handle, 0x102C)  # CL_DEVICE_VENDOR
    version = _cl_device_str(ocl, handle, 0x102F)  # CL_DEVICE_VERSION
    dtype = _cl_device_type(ocl, handle)
    return {
        "name": name[:80],
        "vendor": vendor[:80],
        "version": version[:80],
        "type": dtype,
    }


def _cl_device_str(ocl: Any, handle: int, param: int) -> str:
    fn = ocl.clGetDeviceInfo
    fn.restype = ctypes.c_int32
    fn.argtypes = [
        ctypes.c_void_p,
        ctypes.c_uint32,
        ctypes.c_size_t,
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_size_t),
    ]
    size = ctypes.c_size_t(0)
    err = int(fn(handle, param, 0, None, ctypes.byref(size)))
    if err != 0 or int(size.value) <= 1:
        return ""
    buf = ctypes.create_string_buffer(int(size.value))
    err = int(fn(handle, param, int(size.value), buf, None))
    if err != 0:
        return ""
    return buf.value.decode("utf-8", errors="replace").strip()


def _cl_device_type(ocl: Any, handle: int) -> str:
    fn = ocl.clGetDeviceInfo
    fn.restype = ctypes.c_int32
    fn.argtypes = [
        ctypes.c_void_p,
        ctypes.c_uint32,
        ctypes.c_size_t,
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_size_t),
    ]
    value = ctypes.c_uint64(0)
    size = ctypes.c_size_t(0)
    err = int(fn(handle, 0x1000, ctypes.sizeof(value), ctypes.byref(value), ctypes.byref(size)))
    if err != 0:
        return "unknown"
    bits = int(value.value)
    labels = []
    if bits & 0x2:
        labels.append("CPU")
    if bits & 0x4:
        labels.append("GPU")
    if bits & 0x8:
        labels.append("ACCELERATOR")
    return ",".join(labels) or "other"


def known_cpu_icd_installed() -> Optional[str]:
    """Caminho de intelocl64.dll em localizacoes oficiais Intel, se existir.

    Args:
        Nenhum.

    Returns:
        Caminho absoluto ou None.

    Raises:
        Nenhum.
    """
    program_x86 = os.environ.get("ProgramFiles(x86)") or ""
    program = os.environ.get("ProgramFiles") or ""
    candidates = [
        os.path.join(
            program_x86, "Common Files", "Intel", "OpenCL", "windows", "intelocl64.dll"
        ),
        os.path.join(program_x86, "Intel", "OpenCL", "intelocl64.dll"),
        os.path.join(program, "Intel", "OpenCL", "intelocl64.dll"),
        os.path.join(
            program_x86,
            "Common Files",
            "Intel",
            "shared files",
            "intel64",
            "intelocl64.dll",
        ),
    ]
    for path in candidates:
        if path and os.path.isfile(path):
            return path
    return None


def ocl_icd_filenames_is_hanging(value: str) -> bool:
    """True se OCL_ICD_FILENAMES aponta para um ICD conhecido por pendurar.

    Args:
        value: Conteudo da variavel.

    Returns:
        bool.

    Raises:
        Nenhum.
    """
    text = str(value or "").strip().lower().replace("\\", "/")
    if not text:
        return False
    base = os.path.basename(text)
    return base in HANGING_ICD_BASENAMES


def subprocess_environ(base: Optional[dict] = None) -> dict:
    """Ambiente para Cas-OFFinder: nao injeta ICD de GPU que pendura.

    Args:
        base: os.environ ou copia; default os.environ.

    Returns:
        dict de ambiente.

    Raises:
        Nenhum.
    """
    env = dict(base if base is not None else os.environ)
    current = str(env.get("OCL_ICD_FILENAMES") or "")
    if ocl_icd_filenames_is_hanging(current):
        env.pop("OCL_ICD_FILENAMES", None)
        env["HELIXSCOPE_OPENCL_ICD_STRIPPED"] = "hanging_gpu_icd"
    cpu = known_cpu_icd_installed()
    if cpu and os.path.isfile(cpu) and not env.get("OCL_ICD_FILENAMES"):
        env["OCL_ICD_FILENAMES"] = cpu
        env["HELIXSCOPE_OPENCL_ICD_SOURCE"] = "official_intel_cpu_icd_on_disk"
    return env


def select_cas_offinder_device(platforms: Optional[Mapping[str, Any]] = None) -> dict:
    """Escolhe C se houver CPU OpenCL; senao G se houver GPU. Sem ICD pendurado.

    Args:
        platforms: Saida de enumerate_platforms; se None, enumera agora.

    Returns:
        Dict device (C|G|A), reason, cpu_available, gpu_available.

    Raises:
        Nenhum.
    """
    data = dict(platforms) if platforms is not None else enumerate_platforms()
    cpu = False
    gpu = False
    accelerator = False
    names: list[str] = []
    for plat in list(data.get("platforms") or []):
        for dev in list(plat.get("devices") or []):
            dtype = str(dev.get("type") or "").upper()
            names.append(f"{dev.get('name')} ({dtype})")
            if "CPU" in dtype:
                cpu = True
            if "GPU" in dtype:
                gpu = True
            if "ACCELERATOR" in dtype:
                accelerator = True
    if cpu:
        device = "C"
        reason = "CPU OpenCL device is available. Cas-OFFinder device C is used."
    elif gpu:
        device = "G"
        reason = (
            "No CPU OpenCL ICD is registered. A GPU OpenCL device is available; "
            "Cas-OFFinder device G is used. GPU is not a requirement."
        )
    elif accelerator:
        device = "A"
        reason = "Only an OpenCL accelerator was found. Device A is used."
    else:
        device = "C"
        reason = (
            "No OpenCL devices found. Device C will be requested and Cas-OFFinder "
            "will fail honestly if the ICD is missing."
        )
    return {
        "device": device,
        "reason": reason,
        "cpu_available": cpu,
        "gpu_available": gpu,
        "accelerator_available": accelerator,
        "devices_seen": names[:12],
        "n_platforms": int(data.get("n_platforms") or 0),
    }


def cas_offinder_opencl_ready(platforms: Optional[Mapping[str, Any]] = None) -> bool:
    """True se existir pelo menos um dispositivo OpenCL enumeravel.

    Args:
        platforms: Saida de enumerate_platforms; se None, enumera agora.

    Returns:
        bool.

    Raises:
        Nenhum.
    """
    data = dict(platforms) if platforms is not None else enumerate_platforms()
    if int(data.get("n_platforms") or 0) <= 0:
        return False
    for plat in list(data.get("platforms") or []):
        if list(plat.get("devices") or []):
            return True
    return False


_DIAGNOSE_CACHE: Optional[dict] = None


def diagnose(*, refresh: bool = False) -> dict:
    """Relatorio completo, sem dados irrelevantes de hardware.

    Args:
        refresh: Se True, reenumera plataformas OpenCL.

    Returns:
        Dict hardware, loader, vendors, platforms, blocker, ready.

    Raises:
        Nenhum.
    """
    global _DIAGNOSE_CACHE
    if _DIAGNOSE_CACHE is not None and not refresh:
        return dict(_DIAGNOSE_CACHE)
    hardware = hardware_summary()
    loader = icd_loader_path()
    vendors = khronos_vendor_icds()
    platforms = enumerate_platforms()
    cpu_icd = known_cpu_icd_installed()
    ready = cas_offinder_opencl_ready(platforms)
    selection = select_cas_offinder_device(platforms)
    blocker = ""
    if not loader:
        blocker = "OpenCL ICD loader (OpenCL.dll / libOpenCL) is missing."
    elif not ready:
        blocker = (
            "No OpenCL devices found. The Windows ICD loader is present but "
            "reports zero devices. Intel CPU Runtime for OpenCL is not installed. "
            "HelixScope will not inject an unregistered GPU ICD via "
            "OCL_ICD_FILENAMES (that path hung Cas-OFFinder 2.4.1 on 2026-08-28)."
        )
    payload = {
        "hardware": hardware,
        "loader_present": bool(loader),
        "loader_basename": os.path.basename(loader) if loader else "",
        "khronos_vendor_count": len(vendors),
        "khronos_vendors": vendors,
        "platforms": platforms,
        "device_selection": selection,
        "official_cpu_icd_path_found": bool(cpu_icd),
        "opencl_ready_for_cas_offinder": ready,
        "preferred_device": selection.get("device") or "C",
        "gpu_required": False,
        "blocker": blocker,
        "official_cpu_runtime": OFFICIAL_INTEL_CPU_RUNTIME_NOTE,
        "official_gpu_runtime": OFFICIAL_INTEL_GPU_NOTE,
        "timestamp": provenance.utc_now(),
        "software_version": provenance.HELIXSCOPE_VERSION,
    }
    _DIAGNOSE_CACHE = payload
    return dict(payload)
