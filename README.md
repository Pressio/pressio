
# Overview

Pressio is a collection of repositories providing reduced-order models (ROMs) capabilties.
Specifically, the whole Pressio framework currently includes the following repositories:

* `pressio`: the main C++ repo containing the ROM algorithms and supporting functionalities;

* `pressio-tutorials`: C++ tutorials explaining how to use `pressio` and its functionalities;

* `pressio-builder`: an auxiliary repo with bash helper scripts for configuring/building/installing `pressio`, and `pressio-tutorials`.

## Questions
For questions, find us on Slack: https://pressioteam.slack.com or open an issue.

## License and Citation
Pressio is released with the following [LICENSE](./LICENSE).

Please see the following axXiv paper: https://arxiv.org/abs/2003.07798

## Structure
For a description of `pressio` code structure, see [here](https://github.com/Pressio/pressio/wiki/Structure-of-pressio).

## Building and Installing

### If you only want to use `pressio` from your code
In this case, since `pressio` is header-only, there is **no building process needed**.
You clone the `pressio` repo, and within your code you include the `pressio/packages` to find the `pressio` headers.
However, since `pressio` uses preprocessor directives to selectively enable/disable code for target TPLs, when you build your code you need to have these preprocessor directives defined.
For example, if your code uses Trilinos, to enabled the Trilinos-related code in `pressio` you need to have `PRESSIO_ENABLE_TPL_TRILINOS` defined *before* you include
the `pressio` headers. The list of CMake options to enable can be found [here](./list_of_cmake_optional_vars_to_enable.md).

### If you want to build the unit and regression tests in `pressio`
Sample cmake configure lines can be found [here](https://github.com/Pressio/pressio/wiki/Sample-CMake-configure-lines-for-pressio).

Follow [this](https://github.com/Pressio/pressio/wiki/Serial-build-of-Pressio-with-tests-enabled) for a basic *serial* build that uses only GTest and Eigen and it is done with `pressio-builder` (which automatically builds) Gtest, Eigen for you.

## Sample codes
While we improve the tutorials, please look at the subdirectory `pressio/tests`

## Disclaimer

* Pressio is work-in-progress. At the time of this writing, it is a fairly young project and things are obviously evolving. Several package would benefit from substantial work on testing and documentation, and this is ongoing. However, `pressio` is functional and has been already tested/deployed on large-scale applications.

## Documentation

The m.css documentation is highly dependent on Doxygen. That being said, we start by initializing doxygen on the main directory:

```
>> cd /path/to/pressio
>> doxygen -g 
```

This creates a default Doxyfile in the main directory. In that Doxfile, we need to make sure to specify the output_directory (which can be set to the pressio directory). If this is left blank, the m.css generator will not work. At this point, you can manually add the files and directories you want doxygen to document. This is done by listing those files and directories after the INPUT tag. You can also set the recursion to True so that you don’t have to list all the documents and directories. You can test to make sure doxygen will produce the documentation by running

```
>> doxygen Doxyfile
```

Which should produce a index.html under the `/path/to/pressio/html directory`. If you open that file up in the browser you will see the standard doxgyen webpage format. The next step will be to use the m.css style to improve this format. 

To compile in the m.css format, we will need a Python conf.py file and a Doxyfile-mcss file. The latter is a simple file that is a go between the Python configure file and the Doxyfile. It should look like this:

```
@INCLUDE                = Doxyfile
GENERATE_HTML           = YES
GENERATE_XML            = YES
XML_PROGRAMLISTING      = NO
```

The conf.py file is where you can change the format of the website (add or remove topics from the nav bar or plugins). An example is the following:

```
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
```

See the m.css documentation for a simple overview of how to specify the nav bar mapping. Essentially, if you want to add a new drop down menu, you have to add a tuple to the LINKS_NAVBAR1 list, in the same format displayed in the conf.py file. In order to add pages, they must first be written in markdown. In this example, we put all the pages in the `/path/to/pressio/`pages folder. You must add this folder and/or documents to the Doxyfile INPUT tag. Doxygen will compile these markdown files and add them into the html directory with the md_pages_*.html format, where * is the name of the file in the pages folder. Make sure MARKDOWN_SUPPORT is set to YES as well.  

That’s basically it. Now, we cd into `/path/to/m.css/documentation` and run `./doxygen.py /path/to/pressio/conf.py.` This will run doxygen with the m.css style. Now, if you open up the html/index.html file you should see the new style. 


*This document is in progress, more details soon.*
