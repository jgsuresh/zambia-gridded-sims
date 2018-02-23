# Goal: do a simple 2D calibration where later runs draw from serialized files of best-fit earlier round.
# Round 0: Long burn-in, and serialize files before intervention (requires serialize_year_first_round)
# Round N>0: (need a start year) use serialized files of best fit from previous run.  Start it 5 years before the serialized_year_first_round
