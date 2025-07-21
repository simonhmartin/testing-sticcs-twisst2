import demes
import demesdraw

import matplotlib.pyplot as plt

# Output SVG.
from matplotlib_inline.backend_inline import set_matplotlib_formats
set_matplotlib_formats("svg")


def admix3(t_C=1e4, t_ABC=1e5, t_ABCD=5e5,                     #split times
                    N_A=1e5, N_C=1e5, N_Bstart = 1e4, N_Bend = 1e4, N_ABC=1e5):   #population sizes
    
    #use the deme builder to set up the demographic history
    b = demes.Builder(time_units="generations")
    b.add_deme("ABC", epochs=[{"start_size":N_ABC,  "end_time":t_ABC}])
    b.add_deme("A",     ancestors=["ABC"],  epochs=[{"start_size":N_A}])
    b.add_deme("C",     ancestors=["ABC"],  epochs=[{"start_size":N_C}])
    b.add_deme("B",     ancestors=["A", "C"], proportions=[0.7, 0.3], start_time=t_C, epochs=[{"start_size":N_Bstart, "end_size":N_Bend}])
    
    graph = b.resolve()
    
    return(graph)




def admix_model(t_C=1e4, t_ABC=5e4, t_ABCD=2e5,                     #split times
                    N_A=1e5, N_B=1e5, N_Cstart = 1e4, N_Cend = 1e4, N_D=1e5, N_ABC=1e5, N_ABCD=1e5):   #population sizes
    
    #use the deme builder to set up the demographic history
    b = demes.Builder(time_units="generations")
    b.add_deme("ABCD",                      epochs=[{"start_size":N_ABCD, "end_time":t_ABCD}])
    b.add_deme("ABC", ancestors=["ABCD"], epochs=[{"start_size":N_ABC,  "end_time":t_ABC}])
    b.add_deme("D",     ancestors=["ABCD"],  epochs=[{"start_size":N_D}])
    b.add_deme("A",     ancestors=["ABC"],  epochs=[{"start_size":N_A}])
    b.add_deme("B",     ancestors=["ABC"],  epochs=[{"start_size":N_B}])
    b.add_deme("C",     ancestors=["A", "B"], proportions=[0.999, 0.001], start_time=t_C, epochs=[{"start_size":N_Cstart, "end_size":N_Cend}])
    
    graph = b.resolve()
    
    return(graph)


graph = admix3()
ax1 = demesdraw.tubes(graph, scale_bar=True, max_time=250000)
ax1.figure.savefig("admix3_model.svg", format="svg")

graph = admix_model()
ax1 = demesdraw.tubes(graph, scale_bar=True, max_time=250000)
ax1.figure.savefig("admix4_model.svg", format="svg")

