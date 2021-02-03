
DOXYFILE = 'Doxyfile-mcss'
MAIN_PROJECT_URL: "Pressio"

LINKS_NAVBAR1 = [
  # Get Started
  #("<a href=\"md_pages_get_started.html\">Get Started</a>",
  ("<a href=>Get Started</a>", # empty href so that get started only has children
   [
     ("<a href=\"md_pages_getstarted_build_and_install.html\">Installation and Dependencies</a>",),
     ("<a href=\"md_pages_getstarted_packages.html\">Packages</a>",),
     ("<a href=\"md_pages_getstarted_build_tests_eigen.html\">Build tests with Eigen only</a>",),
     ("<a href=\"md_pages_getstarted_adapting_app.html\">Step-by-step on adapting an app</a>",),
     ("<a href=\"md_pages_getstarted_adapter_continuous_time_api.html\">Continuous-time API</a>",),
     ("<a href=\"md_pages_getstarted_adapter_discrete_time_api.html\">Discrete-time API</a>",)
   ]),

  # tutorials
  ("<a href=>Tutorials</a>", #("<a href=\"md_pages_tutorials.html\">Tutorials</a>",
   [
     ("<a href=\"md_pages_tutorials_tutorial1.html\">Linear Decoder</a>", ),
     ("<a href=\"md_pages_tutorials_tutorial2.html\">Default Galerkin explicit</a>", ),
   ]),

  # # demos
  # ("<a href=>Demos</a>", #("<a href=\"md_pages_examples.html\">Full Examples</a>",
  #  [
  #    ("<a href=\"md_pages_examples_example1.html\">Example1</a>", )
  #  ]),

  # hyper-reduction
  ("<a href=\"md_pages_hyperreduction.html\">Hyper-reduction</a>",
   [
     ("<a href=\"md_pages_hyperreduction_hyperred_how_to_enable.html\">How to enable hyperreduction</a>",),
     ("<a href=\"md_pages_hyperreduction_hyperred_eigen_example.html\">Eigen example</a>",),
     ("<a href=\"md_pages_hyperreduction_hyperred_tpetra_example.html\">Tpetra example</a>",)
   ]),

  # custom ops
  ("<a href=>Custom Ops</a>", #("<a href=\"md_pages_custom_ops.html\">Custom Ops</a>",
   [
     ("<a href=\"md_pages_custom_ops_default_gal_exp.html\">Ops for Galerkin Explicit Time</a>", ),
   ]),

  # # Adapter API
  # ("<a href=\"md_pages_adapter_api.html\">Adapter API</a>",
  #  [
  #    ("<a href=\"md_pages_adapter_apis_adapter_continuous_time_api.html\">Continuous-time API</a>",),
  #    ("<a href=\"md_pages_adapter_apis_adapter_discrete_time_api.html\">Discrete-time API</a>",)
  #  ]),

  ('Classes', 'annotated', []),
  ('Namespaces', 'namespaces', [])
]

LINKS_NAVBAR2 = []
#   ('Classes', 'annotated', []),
#   ('Namespaces', 'namespaces', [])
# ]

PLUGINS = ['m.htmlsanity', 'm.math', 'm.code', 'm.components', 'm.dot', 'm.images']

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
