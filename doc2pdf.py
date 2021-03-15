# maybe some code for preparation of the document...

from subprocess import run

run(
    [
        "pandoc",
        "README.md",
        "-o",
        "README.pdf",
        "--pdf-engine",
        "weasyprint",
        "--css",
        "documentation/default.css",
    ]
)
