
DOXYFILE = 'Doxyfile-mcss'
MAIN_PROJECT_URL: "Pressio"

LINKS_NAVBAR1 = [

  ("<a href=\"md_pages_introduction.html\">Introduction</a>", []),

  ("<a href=\"md_pages_installation.html\">Installation and Dependencies</a>", []),

  # ("<a href=\"md_pages_components.html\">Components</a>",
  ("<a>Documentation</a>",
   [
     ("<a href=\"md_pages_components_mpl.html\">mpl</a>",),
     ("<a href=\"md_pages_components_utils.html\">utils</a>",),
     ("<a href=\"md_pages_components_type_traits.html\">type_traits</a>",),
     ("<a href=\"md_pages_components_expressions.html\">expressions</a>",),
     ("<a href=\"md_pages_components_ops.html\">ops</a>",),
     ("<a href=\"md_pages_components_qr.html\">qr</a>",),
     ("<a href=\"md_pages_components_linsolvers.html\">linear solvers</a>",),
     ("<a href=\"md_pages_components_nonlinsolvers_general.html\">nonlinear solvers: general</a>",),
     ("<a href=\"md_pages_components_nonlinsolvers_nr.html\">nonlinear solvers: Newton-Raphson</a>",),
     ("<a href=\"md_pages_components_nonlinsolvers_gn.html\">nonlinear solvers: Gauss-Newton</a>",),
     ("<a href=\"md_pages_components_nonlinsolvers_lm.html\">nonlinear solvers: Levenberg-Marquardt</a>",),
     ("<a href=\"md_pages_components_ode_advance.html\">ode advancers</a>",),
     ("<a href=\"md_pages_components_ode_steppers_explicit.html\">ode explicit steppers</a>",),
     ("<a href=\"md_pages_components_ode_steppers_implicit.html\">ode implicit steppers</a>",),
     ("<a href=\"md_pages_components_rom_general.html\">rom: general</a>",),
     ("<a href=\"md_pages_components_rom_decoder.html\">rom: decoder</a>",),
     ("<a href=\"md_pages_components_rom_galerkin.html\">rom: Galerkin </a>",),
     ("<a href=\"md_pages_components_rom_lspg.html\">rom: LSPG</a>",),
     ("<a href=\"md_pages_components_rom_wls.html\">rom: WLS</a>",),
     ("<a href=\"md_pages_components_rom_hyperreduction.html\">rom: hyper-reduction</a>",),
   ]),

  ("<a href=\"md_pages_tutorials.html\">Tutorials</a>", []),

  ("<a href=\"https://github.com/Pressio/pressio\">Github Page</a>", []),

  # # core concepts
  # ("<a>Core Concepts</a>",
  #  [
  #    ("<a href=\"md_pages_coreconcepts_adapter_api.html\">Adapter API</a>",),
  #    ("<a href=\"md_pages_coreconcepts_datatypes.html\">Data Types</a>",),
  #    ("<a href=\"md_pages_coreconcepts_densenature.html\">Dense nature of ROMs</a>",),
  #         ("<a href=\"md_pages_coreconcepts_decoder.html\">ROM-to-FOM Decoder</a>",)
  #    # ("<a href=\"md_pages_coreconcepts_adapting_app.html\">Step-by-step on adapting an app</a>",)
  #  ]),

  # #("<a href=\"md_pages_get_started.html\">Get Started</a>",
  # ("<a>Installation and Dependencies</a>", # empty href so that get started only has children
  #  [
  #    ("<a href=\"md_pages_getstarted_build_and_install.html\">Installation and Dependencies</a>",),
  #    ("<a href=\"md_pages_getstarted_packages.html\">Components</a>",),
  #    ("<a href=\"md_pages_getstarted_build_tests_eigen.html\">Build tests with Eigen only</a>",)
  #  ]),

  #("<a href=>Tutorials</a>", ("<a href=\"md_pages_tutorials.html\">Tutorials</a>",
  # #  [
  # #    ("<a href=\"md_pages_tutorials_tutorial1.html\">Linear Decoder</a>", ),
  # #    ("<a href=\"md_pages_tutorials_tutorial2.html\">Default Galerkin explicit</a>", ),
  # #    ("<a href=\"md_pages_tutorials_tutorial3.html\">LSPG ROM of the shallow water equations</a>", )
  # #  ]),

  # # demos
  # ("<a href=>Demos</a>", #("<a href=\"md_pages_examples.html\">Full Examples</a>",
  #  [
  #    ("<a href=\"md_pages_examples_example1.html\">Example1</a>", )
  #  ]),

  # # custom ops
  # ("<a href=>Custom Ops</a>", #("<a href=\"md_pages_custom_ops.html\">Custom Ops</a>",
  #  [
  #    ("<a href=\"md_pages_custom_ops_default_gal_exp.html\">Ops for Galerkin Explicit Time</a>", ),
  #  ]),

  # # Adapter API
  # ("<a href=\"md_pages_adapter_api.html\">Adapter API</a>",
  #  [
  #    ("<a href=\"md_pages_adapter_apis_adapter_continuous_time_api.html\">Continuous-time API</a>",),
  #    ("<a href=\"md_pages_adapter_apis_adapter_discrete_time_api.html\">Discrete-time API</a>",)
  #  ]),

  # ('Classes', 'annotated', []),
  # ('Namespaces', 'namespaces', [])
]

LINKS_NAVBAR2 = []
#   ('Classes', 'annotated', []),
#   ('Namespaces', 'namespaces', [])
# ]

PLUGINS = ['m.htmlsanity', 'm.math', 'm.code', 'm.components', 'm.dot', 'm.images', 'm.table']

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
