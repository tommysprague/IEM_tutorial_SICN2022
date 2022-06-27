# IEM_tutorial_SICN2022

Tutorial code presented at 2022 Kavli Summer Institute in Cognitive Neuroscience on June 28, 2022

For data - please visit https://osf.io/wtmg2/ and download files into a "data" directory at the same level as the root folder of the repository (note: the OSF repo is from the 2019 version of this tutorial, but the data used is the same)

Tutorial files (the github repository) should be within a folder called IEM_tutorial_SICN2022; at that same level, there should be a folder called "data" with subfolders called "fMRI" and "EEG"

For any questions, please contact organizers Tommy Sprague (github.com/tommysprague) 

We'll cover:
1. How to build an encoding model (`IEM_tutorial_walkthrough.mlx`)
2. How to fit that model to data (`IEM_tutorial_walkthrough.mlx`)
3. How to invert the encoding model to reconstruct channel response profiles (`IEM_tutorial_walkthrough.mlx`)
4. How to _decode_ feature values from channel response profiles (`IEM_tutorial_advanced.m`)

There are lots of topics we won't have time to cover in the tutorial session - these are covered in more detail in `IEM_tutorial_advanced_full.m` 

## FILES:
- `IEM_tutorial_walkthrough.mlx` - interactive MATLAB "live notebook" of tutorial on fundamentals of implementing IEM analyses in fMRI or EEG (includes recommended exercises - will not run without completing these exercises) - look at `IEM_tutorial_walkthrough.pptx` for step-by-step narration
- `IEM_tutorial_walkthrough_withAnswers.mlx` - 'complete' version of above, with 'answers' to exercises - will run as downloaded
- `IEM_tutorial_walkthrough_script.m` - same, as .m file (in case using older version of matlab)

- `IEM_tutorial_advanced.m` - ready-to-go set of advanced analysis routines - should run as downloaded (covered in SICN 2022)
- `IEM_tutorial_advanced_full.m` - even more analyses - includes several topics not addressed during SICN 2022 course

We won't cover this explicitly, but for an example of how to simulate data from multiple conditions and compare results across conditions (and across arbitrary transformations to the basis set), see:
- `IEM_sim_simple.m` - copy of code in github.com/JohnSerences/iem_sim which implements fixed encoding models and differences across conditions, as well as invertible linear transforms applied to basis sets (see Gardner & Liu, 2019, eNeuro and Sprague, Boynton & Serences, 2019, eNeuro)