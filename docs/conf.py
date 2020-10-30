
DOXYFILE = 'Doxyfile-mcss'
MAIN_PROJECT_URL: "Pressio"

LINKS_NAVBAR1 = [
  # Get Started
  #("<a href=\"md_pages_get_started.html\">Get Started</a>",
  ("<a href=>Get Started</a>", # empty href so that get started only has children
   [
     ("<a href=\"md_pages_getstarted_packages.html\">Packages</a>",),
     ("<a href=\"md_pages_getstarted_build_and_install.html\">Installation and Dependencies</a>",),
     ("<a href=\"md_pages_getstarted_adapting_app.html\">Step-by-step on adapting an app</a>",)
   ]),

  # tutorials
  ("<a href=>Tutorials</a>", #("<a href=\"md_pages_tutorials.html\">Tutorials</a>",
   [
     ("<a href=\"md_pages_tutorials_tutorial1.html\">Tutorial1</a>", )
   ]),

  # examples
  ("<a href=>Full Examples</a>", #("<a href=\"md_pages_examples.html\">Full Examples</a>",
   [
     ("<a href=\"md_pages_examples_example1.html\">Example1</a>", )
   ]),

  # hyper-reduction
  ("<a href=\"md_pages_hyperreduction.html\">Hyper-reduction</a>",
   [
     ("<a href=\"md_pages_hyperreduction_hyperred_how_to_enable.html\">How to enable hyperreduction</a>",),
     ("<a href=\"md_pages_hyperreduction_hyperred_eigen_example.html\">Eigen example</a>",),
     ("<a href=\"md_pages_hyperreduction_hyperred_tpetra_example.html\">Tpetra example</a>",)
   ]),

  # Adapter API
  ("<a href=\"md_pages_adapter_api.html\">Adapter API</a>",
   [
     ("<a href=\"md_pages_adapter_apis_adapter_galerkin_api.html\">Galerkin ROM</a>",),
     ("<a href=\"md_pages_adapter_apis_adapter_unsteady_lspg_api.html\">Unsteady LSPG ROM</a>",),
     ("<a href=\"md_pages_adapter_apis_adapter_discrete_time_api.html\">Discrete-time API</a>",)
     #("<a href=\"md_pages_adapter_apis_adapter_steady_lspg_api.html\">Steady LSPG ROM</a>",)
   ])
]

LINKS_NAVBAR2 = [
  ('Classes', 'annotated', []),
  ('Namespaces', 'namespaces', [])
]

PLUGINS = ['m.htmlsanity', 'm.math',
           'm.code', 'm.components',
           'm.dot']

SHOW_UNDOCUMENTED = "YES"

FAVICON = 'favicon.ico'


# STYLESHEETS = [
#     'https://fonts.googleapis.com/css?family=Libre+Baskerville:400,400i,700,700i%7CSource+Code+Pro:400,400i,600',
#     '../css/m-light+documentation.compiled.css'
# ]
# THEME_COLOR = '#91cff4'
# FAVICON = 'favicon-light.png'

# STYLESHEETS = [
#     'https://fonts.googleapis.com/css?family=Source+Sans+Pro:400,400i,600,600i%7CSource+Code+Pro:400,400i,600',
#     '../css/m-dark+documentation.compiled.css'
# ]
# THEME_COLOR = '#cb4b16'
# FAVICON = 'favicon-dark.png'
