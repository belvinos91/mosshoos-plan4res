# OMVS cascade (Manantali -> Gouina -> Felou) - Phase-2 scaffold

Known defect to correct FIRST: the station parameter file assigns Manantali's
head and capacity to the BafingMakana GAUGING STATION with zero active storage,
which removes the reservoir from the model. No cascade result is valid before
this is fixed.

Phase-2 deliverable: coupled 3-plant model where Gouina and Felou (run-of-river)
receive Manantali's release + intermediate inflows; joint water value via the
validated SDDP core with a 2-D (storage, inflow-state) Bellman function.
