DOXYFILE = 'Doxyfile-mcss'
MAIN_PROJECT_URL: "Pressio"


LINKS_NAVBAR1 = [
    ("Introduction", 'pages', [
        ("<a href=\"md_pages_about.html\">About</a>", ),
        ("<a href=\"md_pages_introduction.html\">Introduction</a>",)
    ])
]
LINKS_NAVBAR2 = [
    ('Classes', 'annotated', []),
    ('Modules', 'modules', [])
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