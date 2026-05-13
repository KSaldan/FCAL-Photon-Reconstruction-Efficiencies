This repository is directed at describing the steps, along with usable code, in order to extract the photon reconstruction efficiencies of the FCAL.

# Getting Started
Perhaps the first step is to familiarize with our predicessor's analysis and the methods used to exctract the efficiencies.
Jon Zarlings analysis note can be found here https://halldweb.jlab.org/doc-private/DocDB/ShowDocument?docid=4025, describing two methods used to extract the efficiencies. This uses omega -> 3Pi decay events so we will need both reconstructed data and Monte Carlo for this reaction.
## Final State Reconstruction
The full reaction is as follows ɣp ->(π<sup>+</sup>π<sup>-</sup>π<sup>0</sup>)<sub>⍵</sub>p->π<sup>+</sup>π<sup>-</sup> Ɣ (Ɣ)<sub>m</sub>p.
  * ReactionFilter plugin config: "-PReaction1=1_14_7_8_9_14 -PReaction1:Decay1=7_1_m1 -PReaction1:Flags=F1_B2_TO_M7". Feel       free to mess with the flags to reduce file sizes, they will be huge.
## Flattening the Output Trees
The current flattening process is done using DSelector since there are certain branches that are not included in FSROOT and can be found in the following paths in the HD_UTILITIES repository
`hd_utilities/EfficiencyUtilities/NeutralUtilities/DSelectors`
  * Run the "FCAL_skimmer" in this directory over the files
  * Additional skimming for reduced run time/file size can be done running `hd_utilities/EfficiencyUtilities/NeutralUtilities/AnalyzeOutput/FlattenMissShow.cxx`
## Analysis 
The rest of the analysis code is stored in the `hd_utilities/EfficiencyUtilities/NeutralUtilities/AnalyzeOutput` directory
  * MakeHist.cxx will plot the numerator and denominator for each method and relavant histograms.
  * FitOmegaHists will grab the histograms from the previos code and fit the distributions and calculate the efficienies
Once you familiarize yourself with the analysis, play around with the parameters of each step for better optimization.

## Comparing Different Samples
Repeat the steps of the analysis for different data and Monte Carlo samlples to compare the efficiencies.
