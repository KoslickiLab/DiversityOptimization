import subprocess
import tempfile

# No revcomp
with tempfile.NamedTemporaryFile() as temp_file:
    output_file = temp_file.name

    to_run = f"../src/./Form16SSensingMatrix.py -k 4 -i ../data/97_otus_subset.fasta -o {output_file}"
    res = subprocess.run(to_run, shell=True, stdout=subprocess.DEVNULL)
    if res.returncode != 0:
        temp_file.close()
        raise Exception("Test failed to run")

    known_results ="""0	0	1
0	1	6
0	2	5
0	3	1
0	4	4
0	5	4
0	6	8
0	7	8
0	8	3
0	9	9
0	10	11
0	11	6
0	12	5
0	13	6
0	14	4
0	15	3
0	16	2
0	17	9
0	18	4
0	19	5
0	20	3
0	21	3
0	22	6
0	23	3
0	24	5
0	25	8
0	26	16""".split('\n')

    with open(output_file, 'r') as fid:
        for known_result in known_results:
            line = fid.readline().strip()
            print(line)
            print(known_result)
            assert line == known_result

# With revcomp
with tempfile.NamedTemporaryFile() as temp_file:
    output_file = temp_file.name
    to_run = f"../src/./Form16SSensingMatrix.py -k 4 -i ../data/97_otus_subset.fasta -o {output_file} -c"
    res = subprocess.run(to_run, shell=True, stdout=subprocess.DEVNULL)
    if res.returncode != 0:
        raise Exception("Test failed to run")

    num_lines_check = 10
    known_results ="""0	0	1
0	1	8
0	2	5
0	3	1
0	4	5
0	5	8
0	6	15
0	7	12
0	8	5
0	9	10
0	10	15
0	11	7
0	12	5
0	13	8
0	14	5
0	15	6
0	16	7
0	17	14
0	18	6
0	19	8
0	20	9
0	21	11
0	22	18
0	23	7
0	24	11
0	25	20
0	26	23""".split('\n')

    with open(output_file, 'r') as fid:
        for known_result in known_results:
            line = fid.readline().strip()
            assert line == known_result

print("Tests passed successfully!")

