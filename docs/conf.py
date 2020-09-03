
DOXYFILE = 'Doxyfile-mcss'
MAIN_PROJECT_URL: "Pressio"

LINKS_NAVBAR1 = [
  # user guide
  ("<a href=\"md_pages_user_guide.html\">User Guide</a>",
   [
     ("<a href=\"md_pages_build_and_install.html\">Installation</a>",),
     ("<a href=\"md_pages_pressio_app.html\">Interacting with an app</a>",),
     ("<a href=\"md_pages_prep_app.html\">Step-by-step on adapting an app</a>",)
   ]),
  #
  # tutorials
  ("<a href=\"md_pages_tutorials.html\">Tutorials</a>",
   [
     ("<a href=\"md_pages_tutorial1.html\">Tutorial1</a>", )
   ]),
  #
  # examples
  ("<a href=\"md_pages_examples.html\">Full Examples</a>", []),
  #
  # Adapter API
  ("<a href=\"md_pages_adapter_api.html\">Adapter API</a>",
   [
     ("<a href=\"md_pages_adapter_steady_lspg_api.html\">Steady LSPG ROM</a>",),
     ("<a href=\"md_pages_adapter_unsteady_lspg_api.html\">Unsteady LSPG ROM</a>",)
   ]),
  #
  # various
  ("<a href=\"md_pages_various.html\">Various</a>", 
    [
      ("<a href=\"md_pages_license.html\">License</a>",),
      ("<a href=\"md_pages_formulation_lspg.html\">What is LSPG?</a>",)
    ])
]

LINKS_NAVBAR2 = [
  ('Classes', 'annotated', []),
  ('Namespaces', 'namespaces', [])
]

PLUGINS = ['m.htmlsanity', 'm.math', 'm.code', 'm.components', 'm.dot']

SHOW_UNDOCUMENTED = "YES"


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
