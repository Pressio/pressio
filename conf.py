DOXYFILE = 'Doxyfile-mcss'
MAIN_PROJECT_URL: "Pressio"


LINKS_NAVBAR1 = [
    ("Introduction", 'pages', [
        ("<a href=\"md_pages_test1.html\">About</a>", ),
        ("<a href=\"md_pages_test2.html\">Introduction</a>",)
    ])
]
LINKS_NAVBAR2 = [
    ('Classes', 'annotated', []),
    ('Modules', 'modules', [])
]

PLUGINS = ['m.htmlsanity', 'm.math', 'm.code']

SHOW_UNDOCUMENTED = "YES"