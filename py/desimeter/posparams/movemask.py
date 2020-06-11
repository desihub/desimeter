import yaml
from desimeter.bitmask import BitMask

_bitdefs = yaml.safe_load("""
#- position move mask
movemask:
    - [THETA_STUCK,      0, "Theta motor is not moving"]
    - [PHI_STUCK,        1, "Phi motor is not moving"]
    - [FAILED_FIT,       2, "Calibration fit failed"]
    - [INVALID_TABLE,    3, "I/O issue or empty input table"]
    - [INVALID_AFTER_FILTER,   4, "not enough data in table after filtering"]
""")



movemask = BitMask('movemask', _bitdefs)
