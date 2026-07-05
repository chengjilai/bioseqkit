"""Tests for the command-line interface."""

import json

import pytest

from bioseqkit.cli import main


def test_cli_stats(single_fasta, capsys):
    rc = main(["stats", str(single_fasta)])
    assert rc == 0
    out = json.loads(capsys.readouterr().out)
    assert out["n_seqs"] == 1
    assert out["total_length"] == 10


def test_cli_revcomp(single_fasta, capsys):
    rc = main(["revcomp", str(single_fasta)])
    assert rc == 0
    out = capsys.readouterr().out
    assert "GTACGTACGT" in out


def test_cli_translate(single_fasta, capsys):
    rc = main(["translate", str(single_fasta)])
    assert rc == 0
    out = capsys.readouterr().out
    assert "frame+1" in out
    assert "frame-3" in out


def test_cli_kmer(multiline_fasta, capsys):
    rc = main(["kmer", str(multiline_fasta), "-k", "2", "--top", "3"])
    assert rc == 0
    out = capsys.readouterr().out.strip().splitlines()
    assert len(out) <= 3
    assert all("\t" in line for line in out)


def test_cli_index_and_fetch(tmp_path, capsys):
    fa = tmp_path / "g.fa"
    fa.write_text(">chr1\n" + "ACGT" * 10 + "\n")
    assert main(["index", str(fa)]) == 0
    capsys.readouterr()
    assert main(["fetch", str(fa), "chr1:1-4"]) == 0
    out = capsys.readouterr().out
    assert "ACGT" in out


def test_cli_requires_command(capsys):
    with pytest.raises(SystemExit):
        main([])
