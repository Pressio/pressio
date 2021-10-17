
DOXYFILE = 'Doxyfile-mcss'
MAIN_PROJECT_URL: "Pressio"

LINKS_NAVBAR1 = [

  ("<a href=\"md_pages_introduction.html\">Introduction</a>", []),

  ("<a href=\"md_pages_installation.html\">Installation</a>", []),

  ("<a>Documentation/A</a>",
   [
     ("<a href=\"md_pages_components_mpl.html\">mpl</a>",),
     ("<a href=\"md_pages_components_utils.html\">utils</a>",),
     ("<a href=\"md_pages_components_type_traits.html\">type_traits</a>",),
     ("<a href=\"md_pages_components_expressions.html\">expressions</a>",),
     ("<a href=\"md_pages_components_ops.html\">ops</a>",),
     ("<a href=\"md_pages_components_qr.html\">qr</a>",),
     ("<a href=\"md_pages_components_linsolvers.html\">linear solvers</a>",),

     #
     ("<a> nonlinear solvers: </a>",),
     ("<a href=\"md_pages_components_nonlinsolvers_general.html\"> &emsp; - general info</a>",),
     ("<a href=\"md_pages_components_nonlinsolvers_system_api.html\"> &emsp; - problem class API</a>",),
     ("<a href=\"md_pages_components_nonlinsolvers_nr.html\"> &emsp; - Newton-Raphson</a>",),
     ("<a href=\"md_pages_components_nonlinsolvers_gn.html\"> &emsp; - Gauss-Newton</a>",),
     ("<a href=\"md_pages_components_nonlinsolvers_lm.html\"> &emsp; - Levenberg-Marquardt</a>",),
     ("<a href=\"md_pages_components_nonlinsolvers_irls.html\"> &emsp; - irls</a>",),

     #
     ("<a> ode: </a>",),
     ("<a href=\"md_pages_components_ode_advance.html\"> &emsp; - advancers</a>",),
     ("<a href=\"md_pages_components_ode_steppers_explicit.html\"> &emsp; - explicit steppers</a>",),
     ("<a href=\"md_pages_components_ode_steppers_implicit.html\"> &emsp; - implicit steppers</a>",),
   ]),

  ("<a>Documentation/B</a>",
   [
     ("<a> rom: </a>",),
     #("<a href=\"md_pages_components_rom_general.html\"> &emsp; - general info</a>",),
     ("<a href=\"md_pages_components_rom_fom_apis.html\"> &emsp; - FOM adapter API</a>",),
     ("<a href=\"md_pages_components_rom_decoder.html\">  &emsp; - decoder</a>",),

     #
     ("<a href=\"md_pages_components_rom_galerkin.html\"> &emsp; - Galerkin </a>",),
     ("<a href=\"md_pages_components_rom_galerkin_default.html\"> &emsp; &emsp; - default problem </a>",),
     ("<a href=\"md_pages_components_rom_galerkin_hypred.html\"> &emsp; &emsp; - hyper-reduced problem </a>",),
     ("<a href=\"md_pages_components_rom_galerkin_masked.html\"> &emsp; &emsp; - masked problem </a>",),

     #
     ("<a href=\"md_pages_components_rom_lspg_steady.html\"> &emsp; - LSPG: steady </a>",),
     ("<a href=\"md_pages_components_rom_lspg_default_steady.html\"> &emsp; &emsp; - default problem </a>",),
     ("<a href=\"md_pages_components_rom_lspg_hypred_steady.html\"> &emsp; &emsp; - hyper-reduced problem </a>",),
     ("<a href=\"md_pages_components_rom_lspg_masked_steady.html\"> &emsp; &emsp; - masked problem </a>",),

     #
     ("<a href=\"md_pages_components_rom_lspg_unsteady.html\"> &emsp; - LSPG: unsteady </a>",),
     ("<a href=\"md_pages_components_rom_lspg_default.html\"> &emsp; &emsp; - default problem </a>",),
     ("<a href=\"md_pages_components_rom_lspg_hypred.html\"> &emsp; &emsp; - hyper-reduced problem </a>",),
     ("<a href=\"md_pages_components_rom_lspg_masked.html\"> &emsp; &emsp; - masked problem </a>",),

     #
     ("<a href=\"md_pages_components_rom_wls.html\">rom: WLS</a>",),
   ]),

  ("<a href=\"https://pressio.github.io/pressio-tutorials/html/index.html\">Tutorials</a>", []),
  #("<a href=\"md_pages_tutorials.html\">Tutorials</a>", []),

  ("<a href=\"https://github.com/Pressio/pressio\">Github</a>", []),

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
