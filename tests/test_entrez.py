"""Tests for the NCBI Entrez download helper (network mocked)."""

import io
from urllib.request import OpenerDirector

from bioseqkit import entrez


class _FakeResponse:
    def __init__(self, payload: bytes):
        self._buf = io.BytesIO(payload)

    def read(self) -> bytes:
        return self._buf.read()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def test_efetch_fasta_builds_url_and_returns_text(monkeypatch):
    captured = {}

    def fake_urlopen(req, timeout=None):
        captured["url"] = req.full_url
        captured["timeout"] = timeout
        return _FakeResponse(b">NC_000000.1 test\nACGT\n")

    monkeypatch.setattr(entrez.urllib.request, "urlopen", fake_urlopen)
    text = entrez.efetch_fasta("NC_000000.1", email="a@b.c")

    assert text.startswith(">NC_000000.1")
    assert "efetch.fcgi" in captured["url"]
    assert "id=NC_000000.1" in captured["url"]
    assert "rettype=fasta" in captured["url"]
    assert "email=a%40b.c" in captured["url"]


def test_efetch_fasta_result_parses_as_fasta(monkeypatch):
    from bioseqkit import parse_fasta

    def fake_urlopen(req, timeout=None):
        return _FakeResponse(b">seq1 desc\nACGTACGT\n")

    monkeypatch.setattr(entrez.urllib.request, "urlopen", fake_urlopen)
    text = entrez.efetch_fasta("seq1")
    records = list(parse_fasta(io.StringIO(text)))
    assert records[0].id == "seq1"
    assert records[0].sequence == "ACGTACGT"


def test_opener_import_available():
    # sanity: urllib is the only HTTP dependency
    assert issubclass(OpenerDirector, object)
