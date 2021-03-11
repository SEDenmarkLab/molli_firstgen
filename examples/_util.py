import colorama


class ForeColor:

    COLORS = {
        'red': colorama.Fore.RED,
        'green': colorama.Fore.GREEN,
        'blue': colorama.Fore.BLUE,
        'yellow': colorama.Fore.YELLOW,
        'magenta': colorama.Fore.MAGENTA,
        'default': colorama.Fore.WHITE,
    }

    def __init__(self, color: str = 'default'):
        self.color = color.lower()

    def __enter__(self):
        c = self.__class__.COLORS[self.color]
        print(c, end='')
        # return self

    def __exit__(self, *args):
        print(colorama.Style.RESET_ALL, end='')


if __name__ == "__main__":

    with ForeColor('yellow'):
        print("THIS MODULE ONLY PROVIDES CONVENIENCE")

    print("IT IS NOT MEANT TO BE EXECUTED STANDALONE.")
