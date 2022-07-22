# Octave miscellaneous

Miscellaneous <a href="https://www.gnu.org/software/octave/">Octave</a> code. 
It might also work in <a href="http://mathworks.com">Matlab</a>).

## License
All software of my authorship in this page is published under the <a href="https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html">GNU General Purpose License version 2.</a>

## Description 

<ul>

<li><a href="software/trace_rings_DEMO.m">Trace Rings Demonstration</a> </li>

<br>
<li><a href="software/Windowed_Functions_DEMO.m">Windowed Functions Demonstration</a> 
shows windowed functions and their application in procesing audio, kinetic and visual data. It is completed with the 
<a href="software/window_function.m">window_function</a>  and 
<a href="software/multiplot_structure.m">multiplot_structure</a> functions.
Test data:
<a href="software/Mulla_Sanat_On_puolikertosae.wav">audio (musical audio)</a>, 
<a href="software/square_rest_circle_period4s.wii">kinetic (accelerometer)</a>, 
<a href="software/juan_2005_mexico_GREY.jpg">visual (static image)</a>.
</li>

<br>
<li><a href="software/novelty_score_DEMO.m">Novelty Score Demonstration</a> </li>

<br>
<li><a href="software/self_similarity_matrix_DEMO.m">Self-Similarity Matrix Demonstration</a> </li>

<br>
<li>The following two scripts are intended for clustering and classification of timeseries. They use Dynamic Time Warping as a measure of distance and One Nearest Neighbour as a classifier. Both require the <a href="software/dtw.m">dtw</a> function by <a href="http://quanthu.com">Quan Wang</a>.</li>

<ul>
<li><a href="software/unsupervised_timeseries_clustering.m">Unsupervised Timeseries Clustering</a> is a script that groups timeseries according to their similarity.</li>

<li><a href="software/supervised_timeseries_classification.m">Supervised Timeseries Classification</a> is a script that classifies timeseries according to their similarity with given examples.</li>
</ul>

<br>
<li><a href="software/generate_labelled_timeseries.m">Generate Labelled Timeseries</a>
is a script that makes artificial timeseries resembling bodily movement to music. The timeseries are labelled and can be used as testing input of <a href="software/unsupervised_timeseries_clustering.m">Unsupervised Timeseries Clustering</a> or <a href="software/supervised_timeseries_classification.m">Supervised Timeseries Classification</a>.</li>

<br>
<li><a href="software/dgramcluster.m">dgramcluster</a> is a function that finds the biggest clusters in a hierarchical tree with a heuristic of greatest distance between  nodes. 
It can be used by <a href="software/unsupervised_timeseries_clustering.m">Unsupervised Timeseries Clustering</a> but not necessarily since it is embedded in the code.</li>

<br>
<li><a href="software/binary_sequences_similarity_demo.m">Binary Sequences Similarity Demonstration</a> 
is a script that shows different ways of assessing similarity between binary sequences, 
for example musical segmentation boundaries. </li>

<br>
<li><a href="software/binseqsi.m">Binary Sequences Similarity (version 2)</a> 
is a function that computes the similarity between two binary sequences, 
for example sequences of musical segmentation boundaries. 
It can be used by <a href="software/binary_sequences_similarity_demo.m">Binary Sequences Similarity Demonstration</a> 
but not necessarily since it is embedded in the code.</li>

<br>
<li><a href="software/kernel_density_estimation_demo.m">Kernel Density Estimation Demonstration</a> 
is a script that shows how to estimate density using a kernel method, 
for example to represent in one sequence several sequences of musical segmentation boundaries. </li>

<br>
<li><a href="software/physcorr.m">physcorr</a> 
is a function that computes the physical correlation of two input signals, 
for example to assess the similarity between musical segmentation boundaries that have been smoothed by convolving them with a gaussian kernel. 
It can be used by <a href="software/binary_sequences_similarity_demo.m">Binary Sequences Similarity Demonstration</a> 
but not necessarily since it is embedded in the code.</li>

<br>
<li><a href="software/gaussian_function_demo.m">Gaussian Function Demonstration</a> 
is a script that shows the workings of the Gaussian Function (also known as normal distribution curve), 
which can be used for example as a kernel to smooth a sequence of musical segmentation boundaries. </li>

<br>
<li><a href="software/linear_fit_demo.m">Linear Fit Demonstration</a> 
is a script that shows two methods to fit a regression model into two-dimensional data: polynomial fit and 
the normal equation. </li>

<br>
<li><a href="software/logspiral_demo.m">Logarithmic Spiral Demonstration</a> is a script that makes logarithmic spirals in different ways. Just for fun.</li>

</ul>
