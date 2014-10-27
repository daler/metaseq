from metaseq import results_table
import metaseq
import numpy as np

fn = metaseq.example_filename('ex.deseq')
d = results_table.ResultsTable(fn)


def test_dataframe_access():

    # different ways of accessing get the same data in memory
    assert d.id is d.data.id
    assert d['id'] is d.data.id

def test_dataframe_subsetting():
    assert all(d[:10].data == d.data[:10])
    assert all(d.update(d.data[:10]).data == d.data[:10])

def test_copy():
    e = d.copy()
    e.id = 'a'
    assert e.id[0] == 'a'
    assert d.id[0] != 'a'

def smoke_tests():
    #smoke test for repr
    print repr(d)

def test_db():

    # should work
    d.attach_db(None)

    d.attach_db(metaseq.example_filename('dmel-all-r5.33-cleaned.gff.db'))
