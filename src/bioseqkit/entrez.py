"""Optional NCBI Entrez download helper.

Fetches reference sequences from the NCBI E-utilities ``efetch`` endpoint using
only the standard library (``urllib``), so no third-party HTTP client is
required. Network access is required at call time.
"""

from __future__ import annotations

import urllib.parse
import urllib.request

__all__ = ["efetch_fasta"]

_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def efetch_fasta(
    accession: str,
    db: str = "nuccore",
    email: str | None = None,
    api_key: str | None = None,
    timeout: float = 30.0,
) -> str:
    """Download a FASTA record from NCBI by accession and return it as text."""
    params = {
        "db": db,
        "id": accession,
        "rettype": "fasta",
        "retmode": "text",
    }
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key
    url = f"{_EFETCH}?{urllib.parse.urlencode(params)}"
    req = urllib.request.Request(url, headers={"User-Agent": "bioseqkit/0.1"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:  # noqa: S310
        return resp.read().decode("utf-8")
