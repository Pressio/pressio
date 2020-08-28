DOXYFILE = 'Doxyfile-mcss'
MAIN_PROJECT_URL: "Pressio"


LINKS_NAVBAR1 = [
    ("<a href=\"md_pages_user_guide.html\">User Guide</a>", [
        ("<a href=\"md_pages_build_and_install.html\">Build and Installation</a>", ),
        ("<a href=\"md_pages_basic_usage.html\">Basic Usage</a>",),
        ("<a href=\"md_pages_advanced_usage.html\">Advanced Usage</a>",),
        ("<a href=\"md_pages_API.html\">API</a>",)
    ]),
    ("<a href=\"md_pages_tutorials.html\">Tutorials</a>", [
        ("<a href=\"md_pages_quickstart.html\">Quickstart</a>", ),
        ("<a href=\"md_pages_examples.html\">Examples</a>",)
    ]),
	("<a href=\"md_pages_team.html\">Team</a>",[])
]
LINKS_NAVBAR2 = [
    ('Classes', 'annotated', []),
    ('Namespaces', 'namespaces', [])
]

PLUGINS = ['m.htmlsanity', 'm.math', 'm.code']

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