#!/usr/bin/env python3

class Protocol:

    def __init__(self):
        self.steps = []

    def __iadd__(self, step):
        from nonstdlib import MagicFormatter
        self.steps.append(step | MagicFormatter(level=2))
        return self

    def __str__(self):
        from textwrap import indent
        formatted_steps = []

        for i, step in enumerate(self.steps, 1):
            number = "{}. ".format(i)
            padding = ' ' * len(number)
            formatted_steps.append(
                    indent(number + step, ' ' * len(number)).strip())

        return '\n\n'.join(formatted_steps)
