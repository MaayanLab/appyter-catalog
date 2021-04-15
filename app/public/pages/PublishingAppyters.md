# Publishing Appyters

To publish an Appyter you created on the catalog, please submit pull request on GitHub for the Appyters Catalog GitHub repository (<https://github.com/MaayanLab/appyter-catalog>).

The [example appyter](https://github.com/MaayanLab/appyter-catalog/tree/master/appyters/example) may act as a template, it demonstrates some appyter features and includes the fundamental structure necessary around any given appyter to be published in the catalog.

Appyters should be added to the `appyters` subdirectory in their own directory along with the jupyter notebook for that appyter, the `appyter.json` file, and the relevant python, R, or system dependencies to Dockerize the appyter.

All pull requests are validated automatically by a github action to ensure they have the correct format and run with their defaults prior to manual review. It is possible to run this validation locally with `python3 validate/validate_merge.py` from the root of the git branch with your changes.
