from BaseMountRetrieve.basemountretrieve import *
from pathlib import Path

test_dir = Path(__file__).parent / 'tests'

#
# def test_extract_run_name():
#     run_name = extract_run_name(test_dir / 'SampleSheet.csv')
#     assert run_name == 'test_run_name'
#
#
# def test_get_readpair():
#     sample_id = 'SAMPLE-ID-01'
#     fastq_file_list = [
#         Path('/fake/SAMPLE-ID-01_S1_L001_R1_001.fastq.gz'), Path('/fake/SAMPLE-ID-01_S1_L001_R2_001.fastq.gz'),
#         Path('/fake/SAMPLE-ID-02_S1_L001_R1_001.fastq.gz'), Path('/fake/SAMPLE-ID-02_S1_L001_R2_001.fastq.gz'),
#     ]
#     r1, r2 = get_readpair(sample_id=sample_id,
#                           fastq_file_list=fastq_file_list)
#     assert r1 == Path('/fake/SAMPLE-ID-01_S1_L001_R1_001.fastq.gz')
#     assert r2 == Path('/fake/SAMPLE-ID-01_S1_L001_R2_001.fastq.gz')
#
#
# def test_retrieve_fastqgz():
#     fastq_file_list = retrieve_fastqgz(test_dir)
#     assert len(fastq_file_list) > 0
#
#
# def test_retrieve_sampleids():
#     test_sampleids = [
#         Path('/fake/SAMPLE-ID-01_S1_L001_R1_001.fastq.gz'), Path('/fake/SAMPLE-ID-01_S1_L001_R2_001.fastq.gz'),
#         Path('/fake/SAMPLE-ID-02_S1_L001_R1_001.fastq.gz'), Path('/fake/SAMPLE-ID-02_S1_L001_R2_001.fastq.gz'),
#     ]
#     sampleids = retrieve_sampleids(test_sampleids)
#     assert 'SAMPLE-ID-01' in sampleids
#     assert 'SAMPLE-ID-02' in sampleids
#
#
# def test_populate_sample_dictionary():
#     sample_id_list = ['SAMPLE-ID-01', 'SAMPLE-ID-02']
#     fastq_file_list = [
#         Path('/fake/SAMPLE-ID-01_S1_L001_R1_001.fastq.gz'), Path('/fake/SAMPLE-ID-01_S1_L001_R2_001.fastq.gz'),
#         Path('/fake/SAMPLE-ID-02_S1_L001_R1_001.fastq.gz'), Path('/fake/SAMPLE-ID-02_S1_L001_R2_001.fastq.gz'),
#     ]
#     sample_dictionary = populate_sample_dictionary(sample_id_list, fastq_file_list)
#     assert sample_dictionary['SAMPLE-ID-01'] == [Path('/fake/SAMPLE-ID-01_S1_L001_R1_001.fastq.gz'),
#                                                  Path('/fake/SAMPLE-ID-01_S1_L001_R2_001.fastq.gz')]
#     assert sample_dictionary['SAMPLE-ID-02'] == [Path('/fake/SAMPLE-ID-02_S1_L001_R1_001.fastq.gz'),
#                                                  Path('/fake/SAMPLE-ID-02_S1_L001_R2_001.fastq.gz')]
#
#
# def test_get_sample_dictionary():
#     sample_dictionary = get_sample_dictionary(test_dir)
#     assert sample_dictionary['SAMPLE-ID-01'] == [test_dir / 'SAMPLE-ID-01_S1_L001_R1_001.fastq.gz',
#                                                  test_dir / 'SAMPLE-ID-01_S1_L001_R2_001.fastq.gz']
