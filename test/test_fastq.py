import pytest
from biotools.fastq import FastqReader

PATH = "data/test.fastq"

def test_count_and_mean_length():
    with FastqReader(PATH) as r:
        assert r.count() == 5
        assert r.mean_length() == (10 + 20 + 30 + 10 + 10) / 5   # 16.0

def test_generator_yields_tuples():
    with FastqReader(PATH) as r:
        rec = next(iter(r))
    assert rec == ("read1", "ACGTACGTAC", "IIIIIIIIII")

def test_without_with():
    assert FastqReader(PATH).count() == 5

def test_stats():
    with FastqReader(PATH) as r:
        s = r.compute_stats()
    assert s.n_reads == 5
    assert s.length_counts == {10: 3, 20: 1, 30: 1}
    assert s.mean_quality[0] == 40
    assert s.mean_quality[10] == pytest.approx((39 + 39) / 2)
    assert s.base_content.loc[1, "A"] == pytest.approx(60)

def test_bad_record(tmp_path):
    p = tmp_path / "bad.fastq"
    p.write_text("@r\nACGT\n+\nIII\n")
    with pytest.raises(ValueError):
        FastqReader(str(p)).count()

def test_plots(tmp_path):
    with FastqReader(PATH) as r:
        r.plot_all(str(tmp_path))
    assert len(list(tmp_path.glob("*.png"))) == 3

def test_gzip(tmp_path):
    import gzip, shutil
    gz = tmp_path / "t.fastq.gz"
    with open(PATH, "rb") as src, gzip.open(gz, "wb") as dst:
        shutil.copyfileobj(src, dst)
    with FastqReader(str(gz)) as r:
        assert r.count() == 5
        assert r.mean_length() == 16.0