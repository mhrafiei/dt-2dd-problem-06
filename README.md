# dt-2dd-problem-06

```diff
! READ FIRST! 
```

# A. On MARCC, clone the repository:
>* git config --global --unset http.proxy
>* git config --global --unset https.proxy
>* git config --global user.email "mrafiei1@jhu.edu"
>* git config --global user.name "mrafiei1"
>* git clone https://github.com/mhrafiei/dt-2dd-problem-06.git

# B. Edit master_creator_20200504.m for random locations or master_creator_20200511.m for 6 slip systems and 6 partials, and define different values for disl_num in a list (e.g. data repository for disl_num[0], disl_num[1], ...):
>* cd dt-2dd-problem-06/
>* module load matlab
>* matlab -nodisplay -nosplash -nodesktop -r 'master_creator_20200511;'
>* exit

# C. Submit first round of jobs on MARCC:
>* bash ja_sub_1to1.sh

# D. Once completed, submit second round of jobs on MARCC: 
>* bash ja_sub_2to5.sh

# E. For a particular dislocation number, say 10, copy final .mat and .txt data to the data-matlab folder under ml-2dd-problem-06:
>* cp ~/dt-2dd-problem-06/0000000010/results/data_* ~/ml-2dd-problem-06/data-matlab/
# F. run code_data.py (check the ReadMe corresponding to ml-2dd-problem-06)
