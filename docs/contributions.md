# Appyter Catalog Contribution Guide

## Who should contribute?
Anyone who wants to! Please create a pull request to add new content.

## Quick Start
1. Install the appyter cli and its optional dependency cookiecutter.
  ```bash
  pip install appyter cookiecutter
  ```
2. Create an initial appyter in the `appyters` directory using `appyter init`
  ```bash
  cd appyters
  appyter init
  # repond to prompts
  ```
3. Review and update all generated files in the directory.
  - `*.ipynb`: The jupyter notebook is the source code of your appyter, for more information about developing appyters, refer to the appyter docs: <https://github.com/maayanLab/appyter/>. It might also be helpful to review other appyters on the catalog.
  - `README.md`: is used on the appyter landing page of the catalog, rather than setup instructions it should have a longer, more elaborate, description of your appyter for potential users of it.
  - `appyter.json`: has basic metadata about your appyter used for filtering and searching it in the catalog. You should update it with pertinent tags (preferrably, overlapping tags to those seen on the current catalog), and you should also add an image.
    - image: An image for your appyter should be located in a directory `static` under your appyter directory, it should be `png` and precisely `1280x720px`. The filename should be used for the `appyter.json` `image` key.
  - `requirements.txt`: has the python packages one needs to install to use all features of the appyter
4. Validate your appyter. Note that this is run when you make a Pull Request as well.
  ```bash
  # install dependencies to run the tests
  pip install -r compose/requirements.txt -r validate/requirements.txt
  # run the tests on your changes to the repo
  python validate/validate_merge.py -vvv
  ```
5. If all goes well, submit a pull request and it should validate there too, be merged, and deployed in production for the wider community of users to access! Thank you!

## Example
An example appyter exists featuring a complete appyter using most appyter features for your reference, it's under [appyters/example/](../appyters/example).
