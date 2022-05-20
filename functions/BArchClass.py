class BetaArch():
    def __init__(self, arctype, start, sequence, score):
        self.arctype = arctype
        self.start = start
        self.end = start + len(sequence) - 1
        self.sequence = sequence
        self.strandL = int((len(self.sequence) - len(self.arctype)) / 2)
        self.L = len(sequence)
        self.score = score

pass