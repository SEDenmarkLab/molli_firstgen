"""
    Miscellaneous useful utilities.
"""
import colorama
from datetime import datetime


class ForeColor:

    COLORS = {
        "red": colorama.Fore.RED,
        "green": colorama.Fore.GREEN,
        "blue": colorama.Fore.BLUE,
        "yellow": colorama.Fore.YELLOW,
        "magenta": colorama.Fore.MAGENTA,
        "default": colorama.Fore.WHITE,
    }

    def __init__(self, color: str = "default"):
        self.color = color.lower()

    def __enter__(self):
        c = self.__class__.COLORS[self.color]
        print(c, end="")
        # return self

    def __exit__(self, *args):
        print(colorama.Style.RESET_ALL, end="")


class WCTimer:
    """
    Wall clock timer
    """

    def __init__(self, desc: str = ""):
        self.desc = desc

    def __enter__(self):
        self.start = datetime.now()

    def __exit__(self, *args):
        self.finish = datetime.now()
        print(self)

    def __str__(self):
        if hasattr(self, "finish"):
            return f"-----\nCompleted: {self.desc}\n Wall clock time: {self.finish - self.start}\n-----"
        else:
            return f"-----\nIn progress: {self.desc}\n Wall clock time: {datetime.now() - self.start}\n-----"
