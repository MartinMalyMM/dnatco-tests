import os
import subprocess
from pathlib import Path
import gemmi
import pytest
from contextlib import contextmanager
from tempfile import NamedTemporaryFile
from urllib.parse import urlparse, unquote
from urllib.request import urlopen


__author__ = "Martin Maly"
__maintainer__ = "Martin Maly"
__email__ = "martin.maly@mrc-lmb.cam.ac.uk"


# PDB 3A3A - crystal structure of human selenocystine tRNA, contains insertion codes, low resolution
# PDB 5JZQ - centrosymmetric crystal structure of Z-DNA, contains alternate conformations, ultra high resolution
# PDB 9BKD - large structure from SPA cryoEM - human Pdcd4 bound to the 40S small ribosomal subunit 


_RCSB = "https://files.rcsb.org"
_REDO = "https://pdb-redo.eu/db"


def rcsb_mmcif(code: str):
    return f"{_RCSB}/download/{code}.cif"


def redo_cif(code: str):
    return f"{_REDO}/{code}/{code}_final.cif"


def redo_mtz(code: str):
    return f"{_REDO}/{code}/{code}_final.mtz"


@contextmanager
def download(url: str):
    """
    Downloads a file from the given URL and saves it to a temporary file.
    Yields a string path to the temporary file.
    Use in a with statement to ensure the file is deleted afterwards.
    """
    urlName = unquote(os.path.basename(urlparse(url).path))
    with urlopen(url, timeout=30) as response:
        name = response.headers.get_filename() or urlName
        name = name.strip().replace(" ", "_")
        name = "".join(c for c in name if c.isalnum() or c in "-_.")
        with NamedTemporaryFile(suffix=f"_{name}", delete=False) as temp:
            while chunk := response.read(1_000_000):
                temp.write(chunk)
        path = Path(temp.name).resolve()
        try:
            yield str(path)
        finally:
            path.unlink(missing_ok=True)


def check_mmcif_overall_tags(block: gemmi.cif.Block, expected_tags: list[str]) -> bool:
    """Check that all expected tags exist in the mmcif block"""
    for tag_label in expected_tags:
        tag_label_full = f"_ndb_struct_ntc_overall.{tag_label}"
        assert block.find_pair(tag_label_full), f"Missing tag: {tag_label_full}"


def check_mmcif_table_columns(
        block: gemmi.cif.Block,
        table_name: str,
        expected_columns: list[str],
        expected_ntc_steps_table="") -> bool:
    """
    Check that a loop exists and contains all expected columns.
    If expected_ntc_steps_table is provided, also compare the content of the table.
    """
    for loop_label in expected_columns:
        loop_label_full = f"{table_name}{loop_label}"
        assert block.find_loop(loop_label_full), f"Missing loop tag: {loop_label_full}"

        table = block.find(table_name, expected_columns)
        assert table, f"Missing table: {table_name}"

        if expected_ntc_steps_table and table_name == "_ndb_struct_ntc_step.":
            block_expected = gemmi.cif.read_string(expected_ntc_steps_table)[0]
            table_expected = block_expected.find(table_name, expected_columns)
            for row_actual, row_expected in zip(table, table_expected):
                assert list(row_actual) == list(row_expected), f"Mismatch in ntc_steps_table at row: {row_actual}"


ntc_steps_table_3a3a = """data_3A3A_expected
#
loop_
_ndb_struct_ntc_step.id
_ndb_struct_ntc_step.name
_ndb_struct_ntc_step.PDB_model_number
_ndb_struct_ntc_step.label_entity_id_1
_ndb_struct_ntc_step.label_asym_id_1
_ndb_struct_ntc_step.label_seq_id_1
_ndb_struct_ntc_step.label_comp_id_1
_ndb_struct_ntc_step.label_alt_id_1
_ndb_struct_ntc_step.label_entity_id_2
_ndb_struct_ntc_step.label_asym_id_2
_ndb_struct_ntc_step.label_seq_id_2
_ndb_struct_ntc_step.label_comp_id_2
_ndb_struct_ntc_step.label_alt_id_2
_ndb_struct_ntc_step.auth_asym_id_1
_ndb_struct_ntc_step.auth_seq_id_1
_ndb_struct_ntc_step.auth_asym_id_2
_ndb_struct_ntc_step.auth_seq_id_2
_ndb_struct_ntc_step.PDB_ins_code_1
_ndb_struct_ntc_step.PDB_ins_code_2
 1        3a3a_A_G_1_C_2  1  1  A   1  G  .  1  A   2  C  .  A   1  A   2  .  .  
 2        3a3a_A_C_2_C_3  1  1  A   2  C  .  1  A   3  C  .  A   2  A   3  .  .  
 3        3a3a_A_C_3_C_4  1  1  A   3  C  .  1  A   4  C  .  A   3  A   4  .  .  
 4        3a3a_A_C_4_G_5  1  1  A   4  C  .  1  A   5  G  .  A   4  A   5  .  .  
 5      3a3a_A_G_5_G_5.A  1  1  A   5  G  .  1  A   6  G  .  A   5  A   5  .  A  
 6    3a3a_A_G_5.A_A_5.B  1  1  A   6  G  .  1  A   7  A  .  A   5  A   5  A  B  
 7      3a3a_A_A_5.B_U_6  1  1  A   7  A  .  1  A   8  U  .  A   5  A   6  B  .  
 8        3a3a_A_U_6_G_7  1  1  A   8  U  .  1  A   9  G  .  A   6  A   7  .  .  
 9        3a3a_A_G_7_A_8  1  1  A   9  G  .  1  A  10  A  .  A   7  A   8  .  .  
10        3a3a_A_A_8_U_9  1  1  A  10  A  .  1  A  11  U  .  A   8  A   9  .  .  
11       3a3a_A_U_9_C_10  1  1  A  11  U  .  1  A  12  C  .  A   9  A  10  .  .  
12      3a3a_A_C_10_C_11  1  1  A  12  C  .  1  A  13  C  .  A  10  A  11  .  .  
13      3a3a_A_C_11_U_12  1  1  A  13  C  .  1  A  14  U  .  A  11  A  12  .  .  
14      3a3a_A_U_12_C_13  1  1  A  14  U  .  1  A  15  C  .  A  12  A  13  .  .  
15      3a3a_A_C_13_A_14  1  1  A  15  C  .  1  A  16  A  .  A  13  A  14  .  .  
16      3a3a_A_A_14_G_15  1  1  A  16  A  .  1  A  17  G  .  A  14  A  15  .  .  
17      3a3a_A_G_15_U_16  1  1  A  17  G  .  1  A  18  U  .  A  15  A  16  .  .  
18      3a3a_A_U_16_G_18  1  1  A  18  U  .  1  A  19  G  .  A  16  A  18  .  .  
19      3a3a_A_G_18_G_19  1  1  A  19  G  .  1  A  20  G  .  A  18  A  19  .  .  
20      3a3a_A_G_19_U_20  1  1  A  20  G  .  1  A  21  U  .  A  19  A  20  .  .  
21    3a3a_A_U_20_C_20.A  1  1  A  21  U  .  1  A  22  C  .  A  20  A  20  .  A  
22    3a3a_A_C_20.A_U_21  1  1  A  22  C  .  1  A  23  U  .  A  20  A  21  A  .  
23      3a3a_A_U_21_G_22  1  1  A  23  U  .  1  A  24  G  .  A  21  A  22  .  .  
24      3a3a_A_G_22_G_23  1  1  A  24  G  .  1  A  25  G  .  A  22  A  23  .  .  
25      3a3a_A_G_23_G_24  1  1  A  25  G  .  1  A  26  G  .  A  23  A  24  .  .  
26      3a3a_A_G_24_G_25  1  1  A  26  G  .  1  A  27  G  .  A  24  A  25  .  .  
27      3a3a_A_G_25_U_26  1  1  A  27  G  .  1  A  28  U  .  A  25  A  26  .  .  
28      3a3a_A_U_26_G_27  1  1  A  28  U  .  1  A  29  G  .  A  26  A  27  .  .  
29      3a3a_A_G_27_C_28  1  1  A  29  G  .  1  A  30  C  .  A  27  A  28  .  .  
30      3a3a_A_C_28_A_29  1  1  A  30  C  .  1  A  31  A  .  A  28  A  29  .  .  
31      3a3a_A_A_29_G_30  1  1  A  31  A  .  1  A  32  G  .  A  29  A  30  .  .  
32      3a3a_A_G_30_G_31  1  1  A  32  G  .  1  A  33  G  .  A  30  A  31  .  .  
33      3a3a_A_G_31_C_32  1  1  A  33  G  .  1  A  34  C  .  A  31  A  32  .  .  
34      3a3a_A_C_32_U_33  1  1  A  34  C  .  1  A  35  U  .  A  32  A  33  .  .  
35      3a3a_A_U_33_U_34  1  1  A  35  U  .  1  A  36  U  .  A  33  A  34  .  .  
36      3a3a_A_U_34_C_35  1  1  A  36  U  .  1  A  37  C  .  A  34  A  35  .  .  
37      3a3a_A_C_35_A_36  1  1  A  37  C  .  1  A  38  A  .  A  35  A  36  .  .  
38      3a3a_A_A_36_A_37  1  1  A  38  A  .  1  A  39  A  .  A  36  A  37  .  .  
39      3a3a_A_A_37_A_38  1  1  A  39  A  .  1  A  40  A  .  A  37  A  38  .  .  
40      3a3a_A_A_38_C_39  1  1  A  40  A  .  1  A  41  C  .  A  38  A  39  .  .  
41      3a3a_A_C_39_C_40  1  1  A  41  C  .  1  A  42  C  .  A  39  A  40  .  .  
42      3a3a_A_C_40_U_41  1  1  A  42  C  .  1  A  43  U  .  A  40  A  41  .  .  
43      3a3a_A_U_41_G_42  1  1  A  43  U  .  1  A  44  G  .  A  41  A  42  .  .  
44      3a3a_A_G_42_U_43  1  1  A  44  G  .  1  A  45  U  .  A  42  A  43  .  .  
45      3a3a_A_U_43_A_44  1  1  A  45  U  .  1  A  46  A  .  A  43  A  44  .  .  
46      3a3a_A_A_44_G_45  1  1  A  46  A  .  1  A  47  G  .  A  44  A  45  .  .  
47      3a3a_A_G_45_C_46  1  1  A  47  G  .  1  A  48  C  .  A  45  A  46  .  .  
48      3a3a_A_C_46_U_47  1  1  A  48  C  .  1  A  49  U  .  A  46  A  47  .  .  
49    3a3a_A_U_47_G_47.A  1  1  A  49  U  .  1  A  50  G  .  A  47  A  47  .  A  
50  3a3a_A_G_47.A_U_47.B  1  1  A  50  G  .  1  A  51  U  .  A  47  A  47  A  B  
51  3a3a_A_U_47.B_C_47.C  1  1  A  51  U  .  1  A  52  C  .  A  47  A  47  B  C  
52  3a3a_A_C_47.C_U_47.D  1  1  A  52  C  .  1  A  53  U  .  A  47  A  47  C  D  
53  3a3a_A_U_47.D_A_47.E  1  1  A  53  U  .  1  A  54  A  .  A  47  A  47  D  E  
54  3a3a_A_A_47.E_G_47.F  1  1  A  54  A  .  1  A  55  G  .  A  47  A  47  E  F  
55  3a3a_A_G_47.F_C_47.G  1  1  A  55  G  .  1  A  56  C  .  A  47  A  47  F  G  
56  3a3a_A_C_47.G_G_47.H  1  1  A  56  C  .  1  A  57  G  .  A  47  A  47  G  H  
57  3a3a_A_G_47.H_A_47.I  1  1  A  57  G  .  1  A  58  A  .  A  47  A  47  H  I  
58  3a3a_A_A_47.I_C_47.J  1  1  A  58  A  .  1  A  59  C  .  A  47  A  47  I  J  
59  3a3a_A_C_47.J_A_47.K  1  1  A  59  C  .  1  A  60  A  .  A  47  A  47  J  K  
60  3a3a_A_A_47.K_G_47.L  1  1  A  60  A  .  1  A  61  G  .  A  47  A  47  K  L  
61    3a3a_A_G_47.L_A_48  1  1  A  61  G  .  1  A  62  A  .  A  47  A  48  L  .  
62      3a3a_A_A_48_G_50  1  1  A  62  A  .  1  A  63  G  .  A  48  A  50  .  .  
63      3a3a_A_G_50_U_51  1  1  A  63  G  .  1  A  64  U  .  A  50  A  51  .  .  
64      3a3a_A_U_51_G_52  1  1  A  64  U  .  1  A  65  G  .  A  51  A  52  .  .  
65      3a3a_A_G_52_G_53  1  1  A  65  G  .  1  A  66  G  .  A  52  A  53  .  .  
66      3a3a_A_G_53_U_54  1  1  A  66  G  .  1  A  67  U  .  A  53  A  54  .  .  
67      3a3a_A_U_54_U_55  1  1  A  67  U  .  1  A  68  U  .  A  54  A  55  .  .  
68      3a3a_A_U_55_C_56  1  1  A  68  U  .  1  A  69  C  .  A  55  A  56  .  .  
69      3a3a_A_C_56_A_57  1  1  A  69  C  .  1  A  70  A  .  A  56  A  57  .  .  
70      3a3a_A_A_57_A_58  1  1  A  70  A  .  1  A  71  A  .  A  57  A  58  .  .  
71      3a3a_A_A_58_U_59  1  1  A  71  A  .  1  A  72  U  .  A  58  A  59  .  .  
72      3a3a_A_U_59_U_60  1  1  A  72  U  .  1  A  73  U  .  A  59  A  60  .  .  
73      3a3a_A_U_60_C_61  1  1  A  73  U  .  1  A  74  C  .  A  60  A  61  .  .  
74      3a3a_A_C_61_C_62  1  1  A  74  C  .  1  A  75  C  .  A  61  A  62  .  .  
75      3a3a_A_C_62_A_63  1  1  A  75  C  .  1  A  76  A  .  A  62  A  63  .  .  
76      3a3a_A_A_63_C_64  1  1  A  76  A  .  1  A  77  C  .  A  63  A  64  .  .  
77      3a3a_A_C_64_C_66  1  1  A  77  C  .  1  A  78  C  .  A  64  A  66  .  .  
78      3a3a_A_C_66_U_67  1  1  A  78  C  .  1  A  79  U  .  A  66  A  67  .  .  
79    3a3a_A_U_67_U_67.A  1  1  A  79  U  .  1  A  80  U  .  A  67  A  67  .  A  
80  3a3a_A_U_67.A_U_67.B  1  1  A  80  U  .  1  A  81  U  .  A  67  A  67  A  B  
81    3a3a_A_U_67.B_C_68  1  1  A  81  U  .  1  A  82  C  .  A  67  A  68  B  .  
82      3a3a_A_C_68_G_69  1  1  A  82  C  .  1  A  83  G  .  A  68  A  69  .  .  
83      3a3a_A_G_69_G_70  1  1  A  83  G  .  1  A  84  G  .  A  69  A  70  .  .  
84      3a3a_A_G_70_G_71  1  1  A  84  G  .  1  A  85  G  .  A  70  A  71  .  .  
85      3a3a_A_G_71_C_72  1  1  A  85  G  .  1  A  86  C  .  A  71  A  72  .  .  """


ntc_steps_table_5jzq = """data_5JZQ_expected
#
loop_
_ndb_struct_ntc_step.id
_ndb_struct_ntc_step.name
_ndb_struct_ntc_step.PDB_model_number
_ndb_struct_ntc_step.label_entity_id_1
_ndb_struct_ntc_step.label_asym_id_1
_ndb_struct_ntc_step.label_seq_id_1
_ndb_struct_ntc_step.label_comp_id_1
_ndb_struct_ntc_step.label_alt_id_1
_ndb_struct_ntc_step.label_entity_id_2
_ndb_struct_ntc_step.label_asym_id_2
_ndb_struct_ntc_step.label_seq_id_2
_ndb_struct_ntc_step.label_comp_id_2
_ndb_struct_ntc_step.label_alt_id_2
_ndb_struct_ntc_step.auth_asym_id_1
_ndb_struct_ntc_step.auth_seq_id_1
_ndb_struct_ntc_step.auth_asym_id_2
_ndb_struct_ntc_step.auth_seq_id_2
_ndb_struct_ntc_step.PDB_ins_code_1
_ndb_struct_ntc_step.PDB_ins_code_2
 1    5jzq_A_DC.A_1_DG.A_2  1  1  A  1  DC  A  1  A  2  DG  A  A   1  A   2  .  .  
 2    5jzq_A_DC.B_1_DG.B_2  1  1  A  1  DC  B  1  A  2  DG  B  A   1  A   2  .  .  
 3    5jzq_A_DG.A_2_DC.A_3  1  1  A  2  DG  A  1  A  3  DC  A  A   2  A   3  .  .  
 4    5jzq_A_DG.B_2_DC.B_3  1  1  A  2  DG  B  1  A  3  DC  B  A   2  A   3  .  .  
 5    5jzq_A_DC.A_3_DG.A_4  1  1  A  3  DC  A  1  A  4  DG  A  A   3  A   4  .  .  
 6    5jzq_A_DC.B_3_DG.B_4  1  1  A  3  DC  B  1  A  4  DG  B  A   3  A   4  .  .  
 7    5jzq_A_DG.A_4_DC.A_5  1  1  A  4  DG  A  1  A  5  DC  A  A   4  A   5  .  .  
 8    5jzq_A_DG.B_4_DC.B_5  1  1  A  4  DG  B  1  A  5  DC  B  A   4  A   5  .  .  
 9      5jzq_A_DC.A_5_DG_6  1  1  A  5  DC  A  1  A  6  DG  .  A   5  A   6  .  .  
10      5jzq_A_DC.B_5_DG_6  1  1  A  5  DC  B  1  A  6  DG  .  A   5  A   6  .  .  
11    5jzq_B_DC.A_7_DG.A_8  1  1  B  1  DC  A  1  B  2  DG  A  B   7  B   8  .  .  
12    5jzq_B_DC.B_7_DG.B_8  1  1  B  1  DC  B  1  B  2  DG  B  B   7  B   8  .  .  
13    5jzq_B_DG.A_8_DC.A_9  1  1  B  2  DG  A  1  B  3  DC  A  B   8  B   9  .  .  
14    5jzq_B_DG.B_8_DC.B_9  1  1  B  2  DG  B  1  B  3  DC  B  B   8  B   9  .  .  
15   5jzq_B_DC.A_9_DG.A_10  1  1  B  3  DC  A  1  B  4  DG  A  B   9  B  10  .  .  
16   5jzq_B_DC.B_9_DG.B_10  1  1  B  3  DC  B  1  B  4  DG  B  B   9  B  10  .  .  
17  5jzq_B_DG.A_10_DC.A_11  1  1  B  4  DG  A  1  B  5  DC  A  B  10  B  11  .  .  
18  5jzq_B_DG.B_10_DC.B_11  1  1  B  4  DG  B  1  B  5  DC  B  B  10  B  11  .  .  
19  5jzq_B_DC.A_11_DG.A_12  1  1  B  5  DC  A  1  B  6  DG  A  B  11  B  12  .  .  
20  5jzq_B_DC.B_11_DG.B_12  1  1  B  5  DC  B  1  B  6  DG  B  B  11  B  12  .  .  """

@pytest.mark.parametrize(
    "pdb_code,expected_values_minimal,expected_values_precise,expected_ntc_steps_table",
    [
        ("3a3a", {
            "_ndb_struct_ntc_overall.confal_score": 50,
            "_ndb_struct_ntc_overall.confal_percentile": 50,
            "_ndb_struct_ntc_overall.num_classified": 60
            }, {
            "_ndb_struct_ntc_overall.num_steps": 85},
            ntc_steps_table_3a3a
        ),
        ("5jzq", {
            "_ndb_struct_ntc_overall.confal_score": 50,
            "_ndb_struct_ntc_overall.confal_percentile": 50,
            "_ndb_struct_ntc_overall.num_classified": 15
        }, {
            "_ndb_struct_ntc_overall.num_steps": 20},
            ntc_steps_table_5jzq
        ),
        ("9bkd", {
            "_ndb_struct_ntc_overall.confal_score": 40,
            "_ndb_struct_ntc_overall.confal_percentile": 40,
            "_ndb_struct_ntc_overall.num_classified": 1000
        }, {
            "_ndb_struct_ntc_overall.num_steps": 1667},
            ""
        ),
    ],
    ids=['3A3A', '5JZQ', '9BKD']
)
def test_classify_and_write_cif(
        pdb_code: str,
        expected_values_minimal: dict,
        expected_values_precise: dict,
        expected_ntc_steps_table: str
    ):
    executable = "classify_and_write_cif"
    with download(rcsb_mmcif(pdb_code)) as ciffile:
        result = subprocess.run(
            [str(executable), str(ciffile)],
            capture_output=True,
            text=True,
            check=True,
        )

    assert result.returncode == 0
    assert result.stderr == ""
    output = gemmi.cif.read_string(result.stdout)
    assert output

    # Check that the mmcif output contains the expected tags
    block = output[0]
    overall_tags = [
        "entry_id",
        "confal_score",
        "confal_percentile",
        "ntc_version",
        "cana_version",
        "num_steps",
        "num_classified",
        "num_unclassified",
        "num_unclassified_rmsd_close",
    ]
    check_mmcif_overall_tags(block, overall_tags)

    expected_loops = {
        "_ndb_struct_ntc_step.":
        [
            "id",
            "name",
            "PDB_model_number",
            "label_entity_id_1",
            "label_asym_id_1",
            "label_seq_id_1",
            "label_comp_id_1",
            "label_alt_id_1",
            "label_entity_id_2",
            "label_asym_id_2",
            "label_seq_id_2",
            "label_comp_id_2",
            "label_alt_id_2",
            "auth_asym_id_1",
            "auth_seq_id_1",
            "auth_asym_id_2",
            "auth_seq_id_2",
            "PDB_ins_code_1",
            "PDB_ins_code_2",
        ],
        "_ndb_struct_ntc_step_summary.":
        [
            "step_id",
            "assigned_CANA",
            "assigned_NtC",
            "confal_score",
            "euclidean_distance_NtC_ideal",
            "cartesian_rmsd_closest_NtC_representative",
            "closest_CANA",
            "closest_NtC",
            "closest_step_golden",
        ],
        "_ndb_struct_ntc_step_parameters.":
        [
            "step_id",
            "tor_delta_1",
            "tor_epsilon_1",
            "tor_zeta_1",
            "tor_alpha_2",
            "tor_beta_2",
            "tor_gamma_2",
            "tor_delta_2",
            "tor_chi_1",
            "tor_chi_2",
            "dist_NN",
            "dist_CC",
            "tor_NCCN",
            "diff_tor_delta_1",
            "diff_tor_epsilon_1",
            "diff_tor_zeta_1",
            "diff_tor_alpha_2",
            "diff_tor_beta_2",
            "diff_tor_gamma_2",
            "diff_tor_delta_2",
            "diff_tor_chi_1",
            "diff_tor_chi_2",
            "diff_dist_NN",
            "diff_dist_CC",
            "diff_tor_NCCN",
            "confal_tor_delta_1",
            "confal_tor_epsilon_1",
            "confal_tor_zeta_1",
            "confal_tor_alpha_2",
            "confal_tor_beta_2",
            "confal_tor_gamma_2",
            "confal_tor_delta_2",
            "confal_tor_chi_1",
            "confal_tor_chi_2",
            "confal_dist_NN",
            "confal_dist_CC",
            "confal_tor_NCCN",
            "details"
        ],
        "_ndb_struct_sugar_step_parameters.":
        [
            "step_id",
            "P_1",
            "tau_1",
            "Pn_1",
            "P_2",
            "tau_2",
            "Pn_2",
            "nu_1_1",
            "nu_1_2",
            "nu_1_3",
            "nu_1_4",
            "nu_1_5",
            "nu_2_1",
            "nu_2_2",
            "nu_2_3",
            "nu_2_4",
            "nu_2_5",
            "diff_nu_1_1",
            "diff_nu_1_2",
            "diff_nu_1_3",
            "diff_nu_1_4",
            "diff_nu_1_5",
            "diff_nu_2_1",
            "diff_nu_2_2",
            "diff_nu_2_3",
            "diff_nu_2_4",
            "diff_nu_2_5"
        ],
    }
    # Check tables and their columns
    for table, columns in expected_loops.items():
        check_mmcif_table_columns(block, table, columns, expected_ntc_steps_table)

    # Check pairs with expected values
    for tag, value in expected_values_minimal.items():
        assert block.find_pair(tag)[1] >= str(value)
    for tag, value in expected_values_precise.items():
        assert block.find_pair(tag)[1] == str(value)