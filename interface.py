"""
LIPP Kit
Interface

This module specifies the instructions interface which all other parts of the toolkit
assume for how parameter settings ought to be communicated. All this is handled
in one place, the constructor of the Instructions class, to ensure persistent
compatibility of serialisation and deserialisation while this standard evolves.

@author: Mario
"""

class Anything:
    def __contains__(self, other):
        return True

class Instructions:
    min_args = 2

    def __init__(self, from_argv: bool):
        """
        Build an Instructions object, either from command line arguments or
        from user input.

        Parameters
        ----------
        from_argv : bool
            If true, builds an Instructions object with attributes from CLAs.
            This should be used by the solver.
            If false, builds the list from which to reconstruct the former.
            This should be used by the controller.

        Returns
        -------
        None.

        """
        self.slist: list[str] = list()

        def get_next(query: str, options: tuple[str]) -> str:
            out: str
            nonlocal from_argv
            if from_argv:
                nonlocal settings_iter
                out = next(settings_iter)
            else:
                q = query + "\n"
                if type(options) is not Anything:
                    for i in range(len(options)):
                        q += str(i) + " - " + options[i] + "\n"
                out = input(q)
                while (type(options) is not Anything) and (out not in tuple(str(i) for i in range(len(options)))):
                    out = input("Invalid response, choose again: ")
                if type(options) is not Anything:
                    out = options[int(out)]
            assert out in options
            return out
        if from_argv:
            from sys import argv
            assert len(argv) >= Instructions.min_args
            self.graphfile, settings_iter = argv[1], iter(argv[2:])
            from os.path import exists
            assert exists('graphs/' + self.graphfile)
        # TODO: build the actual branching settings configuration, using
        # get_next to obtain values and assign them to attributes and append
        # them to slist, which the controller may extract.
        s = get_next("May the path start anywhere, or at a fixed vertex?", ("anywhere", "at a fixed vertex"))
        self.slist.append(s)
        self.has_start = (s == "at a fixed vertex")
        if self.has_start:
            s = get_next("What's the name of the start vertex?", Anything())
            self.slist.append(s)
            self.start = s
        s = get_next("May the path end anywhere, or at a fixed vertex?", ("anywhere", "at a fixed vertex"))
        self.slist.append(s)
        self.has_goal = (s == "at a fixed vertex")
        if self.has_goal:
            s = get_next("What's the name of the goal vertex?", Anything())
            self.slist.append(s)
            self.goal = s
        pass
