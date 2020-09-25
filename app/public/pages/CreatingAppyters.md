# Creating Appyters

Users can contribute their pipelines to the Appyters Catalog by authoring an Appyter. To author an Appyter, you will need to convert your Jupyter Notebook into an Appyter. This can be done by adding Jupyter Notebook “magics” to your notebook. The Appyter module can enable you to turn your Jupyter Notebook into a jinja2 template-driven web application.

The first step is to install the Appyter package
```bash
# Install package from github repository master
pip3 install --user --upgrade git+git://github.com/Maayanlab/appyter.git
```

The Appyter package enables you to serve that notebook on an executable webapp. Once the package is installed, simply type on the command line:
```bash
appyter jupyter_notebook.ipynb
```

If Appyter is not in your PATH environment variable, you can alternatively type:
```bash
python3 -m appyter jupyter_notebook.ipynb
```

A dotenv file (.env) or environment variables can be use to configure HOST, PORT, and PREFIX of the webserver.

Appyters have several mechanisms for extension. Some of these extensions involve built-in features, and others involve overriding or extending built-in features. “Profiles” are template presets for the default application-defined fields; these enable quick beautification of the Appyter with little effort. “Extras” are feature flags that can be used to enable certain opt-in features such as table of contents, or code toggling. The “extra” independent of the profile. But all existing fields or pages can be overridden or extended by means of a documented directory and file structure; overrides placed in the proper location are automatically loaded by appyters making it quite easy to define new fields or fine-tune the application styling without having to make modifications to appyter itself. Static files, appyter fields, jinja2 filters, jinja2 templates, and even flask blueprints or dash apps can all be defined, integrated, extended, and overridden.

Some pre-configured profiles can be used for styling the form (`--profile=profile_name`) see `appyter/profiles`

In debug mode (`--debug`), changes to the notebook will automatically update the webapp.

- Custom fields can be used by putting them in the directory of execution with the following format:  
  `./fields/YourField.py`: Python Field Implementation
- The templates used natively by the application can be modified to provide your own look and feel:  
  `./templates/head.jinja2`: Custom head (title, CSS, scripts)`  
  `./templates/form.jinja2`: Custom form handling`  
  `./templates/fields/StringField.jinja2`: Override field style
- Custom jinja2 filters can be added if necessary  
  `./filters/your_filter.py`: Python jinja2 filter (function)
- Custom externally-referenced resources (i.e. images) can be put under the static directory  
  `./static/img/your_image.png`: Reference in templates with {% static 'img/your_image.png' %}
- Custom blueprints to additionally mount with flask  
  `./blueprints/your_app/__init__.py`: `your_app = Blueprint('your_app', __name__)`

### Modifying the notebook to become an Appyter

The Appyter “magics” can be used to directly serialize and subsequently execute jinja2-style template syntax. This syntax permits a wide range of branches which enable the notebook’s code to be adjusted as needed based on declared “Fields”. These fields represent the type of input field to be used to construct that template variable. Fields are available automatically for all major data types and can be extended to support more specific use cases of input form components. These fields can be extracted by inspection and are eventually used for the purpose of rendering a web form that is the initial user interface of the Appyter.

Create a standard python Jupyter Notebook, make the first cell:
```python
#%%appyter init
from appyter import magic
magic.init(lambda _=globals: _())
```

Normal cells are allowed, you also have access to jinja2-rendered cells:
```python
%%appyter {cell_type}
# ...
```

Supported cell_types:

- `markdown`: Substitute jinja2 template variables, render as a markdown cell
- `hide`: Substitute jinja2 template variables, show it rendered in your notebook with the default values, but when executing publicly, do not show/execute the cell.
- `code_exec`: Substitute jinja2 template variables, render as python, show it rendered in your notebook with the default values and execute it
- `code_eval`: Substitute jinja2 template variables, render as python, show it rendered in your notebook with the default values and execute it, “eval” the last line of the cell and display the result.

By submitting the form for execution, the Appyter assembles all the necessary variables to fully serialize a customized instance of the target instantiation of the template Jupyter Notebook. A file input form field also exists. It facilitates uploading user-submitted files to utilize for a given analysis.
The appyter command line interface (CLI) can be used to easily interact with the appyter feature-set including locating and describing available fields, profiles, and extras as well as facilitating the inspection, construction, evaluation, and serving (via flask) of appyters. Furthermore, the interface facilitates interacting with remote appyter instances using both the appyter REST API and with websockets, allowing inspection and real-time asynchronous evaluation of public appyter endpoints straight from the command line or indeed, as part of a workflow.
