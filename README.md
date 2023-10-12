# pps_project
LSAMP project describing PPS using computational methods

1. MD sim of PPS morphology will be predicted using cmelab/jankflow to write
code based on hoomd to run a simulation.
2. Characterizing PPS by running several different amounts of molecules (N).
A plot will be created of TPS (timesteps per second) vs N.
TPS tells us how fast our simulation is progressing per second.
Number of PPS molecules will exponentially increase at a rate of e^10,
which will require help from Fry, a research computing cluster. This plot will
be written using matplotlib.
3. Run cmeutils function located in cmelab/cmelabutils/structure.py called
gsdRDF() and pass in a gsdfile of PPS simulated, bin_centers vs  RDF.RDF.
RDF curves tell us how likely at r we are to find a neighbor at that
distance.
