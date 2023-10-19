# // geophysical_notes //

My collection of geophysical notes written as Jupyter notebooks.

(Mostly obsolete now, don't have time for a proper update; for example I stopped wasting my time with jupyter widgets and Panel for interactive apps and I am now using Streamlit).

# seismic petrophysics

Do magic things with well log data.

* [Seismic Petrophysics](http://nbviewer.ipython.org/github/aadm/geophysical_notes/blob/master/seismic_petrophysics.ipynb)
* [Seismic Petrophysics / interactive](http://nbviewer.ipython.org/github/aadm/geophysical_notes/blob/master/seismic_petrophysics_interactive.ipynb)
* [Rock physics modeling and templates](http://nbviewer.ipython.org/github/aadm/geophysical_notes/blob/master/rock_physics_modeling.ipynb)


#### support data for "seismic petrophysics"

* `well_2*.txt`: raw log data from Well 2 of [Quantitative Seismic Interpretation (QSI)](https://srb.stanford.edu/quantitative-seismic-interpretation)
* `qsiwell2.csv`: assembled all the logs from various files
* `qsiwell2_frm.csv`: qsiwell2 + fluid replaced elastic logs
* `qsiwell2_augmented.csv`: barebones well data, only Ip, Vp/Vs and LFC (litho-fluid class log)
* `qsiwell2_synthetic.csv`: synthetic data generated through Monte Carlo simulation, same logs as in `qsiwell2_augmented.csv` (Ip, Vp/Vs and LFC)
* `qsiwell2_dataprep.py`: Python script to assemble all the original QSI files


# seismic and rock physics stuff

How to load and display SEG-Y files, plus some simple ways to play with the data, e.g. extracting amplitude informations, adding noise & filtering. Also, a notebook entirely dedicated to wedge modeling and how to reproduce a couple of figures from scientific publications.

* [Playing with seismic](http://nbviewer.ipython.org/github/aadm/geophysical_notes/blob/master/playing_with_seismic.ipynb)
* [Playing with seismic (interactive)](http://nbviewer.ipython.org/github/aadm/geophysical_notes/blob/master/playing_with_seismic_interactive.ipynb)
* [Amplitude extraction](http://nbviewer.ipython.org/github/aadm/geophysical_notes/blob/master/seismic_amplitude_extraction.ipynb)
* [Wedge modeling for variable angles of incidence](http://nbviewer.ipython.org/github/aadm/geophysical_notes/blob/master/wedge_modeling.ipynb)
* [Notes on spectral decomposition](http://nbviewer.ipython.org/github/aadm/geophysical_notes/blob/master/notes_spec_dec.ipynb)
* [Top Heimdal map, or how to reproduce figure 1 from Avseth et al., 2001](http://nbviewer.ipython.org/github/aadm/geophysical_notes/blob/master/top_heimdal_map.ipynb)
* [AVO projections](http://nbviewer.ipython.org/github/aadm/geophysical_notes/blob/master/avo_projections.ipynb)
* [How to calculate AVO attributes](http://nbviewer.ipython.org/github/aadm/geophysical_notes/blob/master/avo_attributes.ipynb)
* [Elastic Impedance](http://nbviewer.ipython.org/github/aadm/geophysical_notes/blob/master/elastic_impedance.ipynb)
* ["The relationship between reflectivity and elastic impedance", or how to reproduce figure 5.62 from Seismic Amplitude by Simm & Bacon (2014)](http://nbviewer.ipython.org/github/aadm/geophysical_notes/blob/master/relationship-reflectivity-elastic-impedance_Simm-Bacon.ipynb)
* [Notes on anisotropic AVO equations](http://nbviewer.ipython.org/github/aadm/geophysical_notes/blob/master/anisotropic_avo.ipynb)
* [AVO Explorer v2: Interactive AVO and AVO classes explorer](http://nbviewer.ipython.org/github/aadm/geophysical_notes/blob/master/avo_explorer_v2.ipynb): meant to be downloaded and run locally.
* [Simple porosity modeling](http://nbviewer.ipython.org/github/aadm/geophysical_notes/blob/master/simple_porosity_modeling.ipynb): how to model porosity variations and its effects on elastic properties using the concept of pore stiffness invariance.


#### support data for "seismic stuff"

* `16_81_PT1_PR.SGY`, `16_81_PT2_PR.SGY`, `16_81_PT3_PR.SGY`, `31_81_PR.SGY`: 2D lines in SEGY format from the [USGS Alaska dataset](http://energy.usgs.gov/GeochemistryGeophysics/SeismicDataProcessingInterpretation/NPRASeismicDataArchive.aspx)
* `3d_farstack.sgy`, `3d_nearstack.sgy`: 3D cubes from the QSI dataset (see above)
* `Top_Heimdal_subset.txt`: interpreted horizon for the QSI near and far angle cubes
* `ST10010ZC11_FAR_T_CROP.seisnc`, `ST10010ZC11_NEAR_MID_T_CROP.seisnc`, `ST10010ZC11_NEAR_T_CROP.seisnc`: angle stacks from [Equinor's Volve dataset](https://data.equinor.com/dataset/Volve)
* Interpreted horizons for the Volve seismic are in `Horizons_TWT/Official_Horizons`
* `volve10r12-full-twt-sub3d.sgy`, `hor_twt_hugin_fm_top.dat`: small segy and interpreted horizon from [segysak documentation examples](https://github.com/trhallam/segysak) (also from Volve dataset)


# miscellaneous

Other notebook of interest, maybe only tangentially related to geophysics, such as a notebook showing a comparison between colormaps (the dreadful _jet_ against a bunch of better alternatives) and another that uses the well known Gardner's equation as an excuse to practice data fitting in Python.

* [Color palettes](http://nbviewer.ipython.org/github/aadm/geophysical_notes/blob/master/colormaps.ipynb)
* [Inverse Gardner](http://nbviewer.ipython.org/github/aadm/geophysical_notes/blob/master/inverse_gardner.ipynb)



## notes on running python

I used to recommend either [Enthought's Canopy Express]((https://www.enthought.com/products/canopy/)) or [Anaconda](https://www.continuum.io/why-anaconda). I haven't been using Canopy for a while now and I'm very particular about installing stuff on my computers, so right now what I use (and suggest everyone else to do) is to install a subset of Anaconda called [miniconda](http://conda.pydata.org/miniconda.html). Starting from this, it's easy to install the packages you need for your work and nothing else. For example, this is how I setup my system for work:

```
$ conda install numpy scipy pandas matplotlib jupyter scikit-learn scikit-image xarray dask netCDF4 bottleneck
$ conda install -c bokeh colorcet
$ conda install -c conda-forge jupyterlab
```

Then I install some additional packages with `pip`:

```
$ pip install bruges lasio segyio
```

Instead of [integrated environments](https://en.wikipedia.org/wiki/Integrated_development_environment) (for example, [Spyder](https://github.com/spyder-ide/spyder)) I simply use a modern (and free!) editor like [VSCode](https://code.visualstudio.com/) ([Atom](https://atom.io/) is a good alternative) to code and write. However, [JupyterLab](https://github.com/jupyterlab/jupyterlab) gets better everyday and it can already be used to do everything in a browser window (but to me it's still slower than a text editor and a jupyter console window); I like the idea of Juyter Notebooks to distribute commented code and simply as a working tool to make code interact with explanatory text and plots.

### using SEG-Y data

To read and write SEG-Y data in Python you need additional libraries like  [ObsPy](http://obspy.org), [Segpy](https://github.com/sixty-north/segpy) or Equinor's [segyio](https://github.com/equinor/segyio). 

My current favourite is however [segysak](https://github.com/trhallam/segysak), based on segyio. In addition to the qualities of the underlying segyio (a 340 Mb file is read in 1 second, while obspy needs 8), it also adds some very cool features to map horizons to seismic cubes (amplitude extraction!).

```
# timeit results using segyio:
1.11 s ± 17.7 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)

# timeit results using obspy:
8.85 s ± 1.07 s per loop (mean ± std. dev. of 7 runs, 1 loop each)
```


## license

<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">Creative Commons Attribution-NonCommercial 4.0 International License</a>.
