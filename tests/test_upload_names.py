"""Upload names are not paths, and empty uploads are not sequences."""

from __future__ import annotations

import pytest

import app


class _Upload:
    def __init__(self, name: str, data: bytes) -> None:
        self.name = name
        self._data = data

    def getvalue(self) -> bytes:
        return self._data


def test_simple_fasta_name_is_read() -> None:
    text = app._read_uploaded_text(_Upload("sample.fasta", b">s\nACGT\n"))
    assert "ACGT" in text


@pytest.mark.parametrize(
    "name",
    [
        "../../etc/passwd.fasta",
        "..\\..\\etc\\passwd.fasta",
        "C:\\Users\\test.fasta",
        "/tmp/test.fasta",
        "....//....//test.fasta",
        "",
    ],
)
def test_path_like_upload_names_are_rejected(name: str) -> None:
    with pytest.raises(ValueError, match="simple file name"):
        app._read_uploaded_text(_Upload(name, b">s\nACGT\n"))


def test_empty_upload_is_rejected() -> None:
    with pytest.raises(ValueError, match="empty"):
        app._read_uploaded_text(_Upload("empty.fasta", b""))
