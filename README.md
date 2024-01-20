# Octave miscellaneous

<img src="https://github.com/juigmend/octave-miscellaneous/raw/main/juan_2005_mexico_GREY.jpg" width="1200" height="240">

Miscellaneous <a href="https://www.gnu.org/software/octave/">Octave</a> code. 
It might also work in <a href="http://mathworks.com">Matlab</a>.

## Description 

<ul>

<li><A HREF="https://github.com/juigmend/octave-miscellaneous/blob/main/matrices_multiplication_DEMO.m">Multiplication of Matrices</A> 
 is a starter's tutorial, including definition of scalars, matrices and vectors; indexing and basic arithmetic. </li>
<li><A HREF="https://github.com/juigmend/octave-miscellaneous/blob/main/loops_DEMO.m">Loops</A> 
 includes definition and examples of basic logic operators; "while" and "for" loops. </li>
 <li><A HREF="https://github.com/juigmend/octave-miscellaneous/blob/main/find_DEMO.m">Find </A> 
 shows differences between matrices and cells; how to find values in an array and a practical application. </li>
 <li><A HREF="https://github.com/juigmend/octave-miscellaneous/blob/main/conditional_statements_DEMO.m">Conditional Statements </A> 
 shows operators 'if', 'else' and 'ifelse'. </li>
<li><A HREF="https://github.com/juigmend/octave-miscellaneous/blob/main/MIR_win_convo_autocorr_DEMO.m">Music Information Retrieval - Windowing, Convolution and Autocorrelation</A> 
 includes definition and examples applied to Music Information Retrieval, of windowed processes, convolution, correlation, cross-correlation and autocorrelation; extraction of harmonics by FFT,
distance and self-similarity matrices; peak detection and segmentation by novelty. </li>  
<li><A HREF="https://github.com/juigmend/octave-miscellaneous/blob/main/MIR_query_melody_ASS.m">MIR Query-by-Melody Assignment</A> 
 includes challenges to apply the knowledge contained in the files above, into a system that exracts discrete frequency segments from audio files, 
 represents those segments as a sequence of musical notes and compares the sequence with a melody in a database. </li>
 
<br>

<li><a href="https://github.com/juigmend/octave-miscellaneous/blob/main/trace_rings_DEMO.m">Trace Rings Demonstration</a> </li>

<br>
<li><a href="https://github.com/juigmend/octave-miscellaneous/blob/main/Windowed_Functions_DEMO.m">Windowed Functions Demonstration</a> 
shows windowed functions and their application in procesing audio, kinetic and visual data. It is completed with the 
<a href="https://github.com/juigmend/octave-miscellaneous/blob/main/window_function.m">window_function</a>  and 
<a href="https://github.com/juigmend/octave-miscellaneous/blob/main/multiplot_structure.m">multiplot_structure</a> functions.
Test data:
<a href="https://github.com/juigmend/octave-miscellaneous/blob/main/Mulla_Sanat_On_puolikertosae.wav">audio (musical audio)</a>, 
<a href="https://github.com/juigmend/octave-miscellaneous/blob/main/square_rest_circle_period4s.wii">kinetic (accelerometer)</a>, 
<a href="https://github.com/juigmend/octave-miscellaneous/blob/main/juan_2005_mexico_GREY.jpg">visual (static image)</a>.
</li>

<br>
<li><a href="https://github.com/juigmend/octave-miscellaneous/blob/main/novelty_score_DEMO.m">Novelty Score Demonstration</a> </li>

<br>
<li><a href="https://github.com/juigmend/octave-miscellaneous/blob/main/self_similarity_matrix_DEMO.m">Self-Similarity Matrix Demonstration</a> </li>

<br>
<li>The following two scripts are intended for clustering and classification of timeseries. They use Dynamic Time Warping as a measure of distance and One Nearest Neighbour as a classifier. Both require the <a href="https://github.com/juigmend/octave-miscellaneous/blob/main/dtw.m">dtw</a> function by <a href="http://quanthu.com">Quan Wang</a>.</li>

<ul><blockquote>
<li><a href="https://github.com/juigmend/octave-miscellaneous/blob/main/unsupervised_timeseries_clustering.m">Unsupervised Timeseries Clustering</a> is a script that groups timeseries according to their similarity.</li>

<li><a href="https://github.com/juigmend/octave-miscellaneous/blob/main/supervised_timeseries_classification.m">Supervised Timeseries Classification</a> is a script that classifies timeseries according to their similarity with given examples.</li>
</ul>

<br>
<li><a href="https://github.com/juigmend/octave-miscellaneous/blob/main/generate_labelled_timeseries.m">Generate Labelled Timeseries</a>
is a script that makes artificial timeseries resembling bodily movement to music. The timeseries are labelled and can be used as testing input of <a href="https://github.com/juigmend/octave-miscellaneous/blob/main/unsupervised_timeseries_clustering.m">Unsupervised Timeseries Clustering</a> or <a href="https://github.com/juigmend/octave-miscellaneous/blob/main/supervised_timeseries_classification.m">Supervised Timeseries Classification</a>.</li>

<br>
<li><a href="https://github.com/juigmend/octave-miscellaneous/blob/main/dgramcluster.m">dgramcluster</a> is a function that finds the biggest clusters in a hierarchical tree with a heuristic of greatest distance between  nodes. 
It can be used by <a href="https://github.com/juigmend/octave-miscellaneous/blob/main/unsupervised_timeseries_clustering.m">Unsupervised Timeseries Clustering</a> but not necessarily since it is embedded in the code.</li>

<br>
<li><a href="https://github.com/juigmend/octave-miscellaneous/blob/main/binary_sequences_similarity_demo.m">Binary Sequences Similarity Demonstration</a> 
is a script that shows different ways of assessing similarity between binary sequences, 
for example musical segmentation boundaries. </li>

<br>
<li><a href="https://github.com/juigmend/octave-miscellaneous/blob/main/binseqsi.m">Binary Sequences Similarity (version 2)</a> 
is a function that computes the similarity between two binary sequences, 
for example sequences of musical segmentation boundaries. 
It can be used by <a href="https://github.com/juigmend/octave-miscellaneous/blob/main/binary_sequences_similarity_demo.m">Binary Sequences Similarity Demonstration</a> 
but not necessarily since it is embedded in the code.</li>

<br>
<li><a href="https://github.com/juigmend/octave-miscellaneous/blob/main/kernel_density_estimation_demo.m">Kernel Density Estimation Demonstration</a> 
is a script that shows how to estimate density using a kernel method, 
for example to represent in one sequence several sequences of musical segmentation boundaries. </li>

<br>
<li><a href="https://github.com/juigmend/octave-miscellaneous/blob/main/physcorr.m">physcorr</a> 
is a function that computes the physical correlation of two input signals, 
for example to assess the similarity between musical segmentation boundaries that have been smoothed by convolving them with a gaussian kernel. 
It can be used by <a href="https://github.com/juigmend/octave-miscellaneous/blob/main/binary_sequences_similarity_demo.m">Binary Sequences Similarity Demonstration</a> 
but not necessarily since it is embedded in the code.</li>

<br>
<li><a href="https://github.com/juigmend/octave-miscellaneous/blob/main/gaussian_function_demo.m">Gaussian Function Demonstration</a> 
is a script that shows the workings of the Gaussian Function (also known as normal distribution curve), 
which can be used for example as a kernel to smooth a sequence of musical segmentation boundaries. </li>

<br>
<li><a href="https://github.com/juigmend/octave-miscellaneous/blob/main/linear_fit_demo.m">Linear Fit Demonstration</a> 
is a script that shows two methods to fit a regression model into two-dimensional data: polynomial fit and 
the normal equation. </li>

<br>
<li><a href="https://github.com/juigmend/octave-miscellaneous/blob/main/logspiral_demo.m">Logarithmic Spiral Demonstration</a> is a script that makes logarithmic spirals in different ways. Just for fun.</li>

</ul>

## License
All software of my authorship in this page is published under the <a href="https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html">GNU General Purpose License version 2.</a>
