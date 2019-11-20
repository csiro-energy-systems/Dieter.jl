Dieter.jl Change Log
=========================

### 2019-Nov-20 Commit:
- Added logic in code to allow model components (EV, Heat, H2 etc.)
to be switched in or out of model by setting corresponding setting to `missing`
- Corrected deprecated DataFrames usage: `df[!,:colname]` replaces `df[:colname]`
- Added folder path variable to save and merge results data, while defaulting the directory
name to `"single_results"`
- Added function `initialise_data_file_dict!` to initialise a basic file structure
 for accessing the model data.
- Renamed and added this CHANGELOG.
- Added a file (`run_script.jl`) demonstrating a full workflow of the code.


### 2019-Nov-15 Commit:
- Modified structure of data-parsing functions to directly take a
prepared DataFrame rather than have each function parse CSV files.
- Existing code has been updated for latest DataFrame syntax.
- Utility functions have been modified to use this workflow.
