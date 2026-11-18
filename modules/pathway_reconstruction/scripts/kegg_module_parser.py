#!/usr/bin/env python3
"""
kegg_module_parser.py
Recursive AST parser and evaluator for KEGG Module definition grammar.
Supports:
  - Space ' ': AND (sequential required steps)
  - Comma ',': OR (alternative enzymes or parallel paths)
  - Plus '+': Complex (multi-protein complexes where all subunits are required)
  - Minus '-': Optional components (missing does not penalize parent step)
  - Parentheses '()': Arbitrary nested grouping
"""

import re

class Node:
    def evaluate(self, present_kos):
        raise NotImplementedError

class KONode(Node):
    def __init__(self, ko_id):
        self.ko_id = ko_id.upper()

    def evaluate(self, present_kos):
        if self.ko_id in present_kos:
            return 1.0, True, {self.ko_id}, set()
        return 0.0, False, set(), {self.ko_id}

    def __repr__(self):
        return self.ko_id

class OptionalNode(Node):
    def __init__(self, child):
        self.child = child

    def evaluate(self, present_kos):
        c_score, c_done, matched, missing = self.child.evaluate(present_kos)
        # Optional: if present, contributes full; if missing, treated as complete without penalty
        return 1.0, True, matched, set()

    def __repr__(self):
        return f"-{self.child}"

class ComplexNode(Node):
    def __init__(self, children):
        self.children = children # required + optional subunits

    def evaluate(self, present_kos):
        if not self.children:
            return 1.0, True, set(), set()

        required = [c for c in self.children if not isinstance(c, OptionalNode)]
        all_matched = set()
        for c in self.children:
            _, _, matched, _ = c.evaluate(present_kos)
            all_matched.update(matched)

        if required:
            req_evals = [c.evaluate(present_kos) for c in required]
            all_done = all(done for _, done, _, _ in req_evals)
            all_missing = set()
            for _, _, _, missing in req_evals:
                all_missing.update(missing)
            completeness = sum(sc for sc, _, _, _ in req_evals) / len(required)
        else:
            completeness = 1.0
            all_done = True
            all_missing = set()

        return completeness, all_done, all_matched, all_missing

    def __repr__(self):
        return "+".join(str(c) for c in self.children)

class OrNode(Node):
    def __init__(self, children):
        self.children = children

    def evaluate(self, present_kos):
        if not self.children:
            return 0.0, False, set(), set()

        evals = [c.evaluate(present_kos) for c in self.children]

        # 1. Optimal branch: any completed branch satisfies the step
        for sc, done, matched, missing in evals:
            if done:
                return sc, done, matched, set()

        # 2. Partially matched branch (sc > 0): pick the branch closest to completion
        positive_evals = [e for e in evals if e[0] > 0.0]
        if positive_evals:
            positive_evals.sort(key=lambda e: (e[0], -len(e[3])), reverse=True)
            return positive_evals[0]

        # 3. No branch has any matches (all sc == 0.0):
        # Format alternative missing paths with | so all functional pathways are visible
        alt_missing_strs = []
        for _, _, _, missing in evals:
            if missing:
                alt_missing_strs.append("+".join(sorted(missing)))

        if alt_missing_strs:
            combined_missing = {"|".join(alt_missing_strs)}
        else:
            combined_missing = set()

        return 0.0, False, set(), combined_missing

    def __repr__(self):
        return "(" + ",".join(str(c) for c in self.children) + ")"

class AndNode(Node):
    def __init__(self, children):
        self.children = children

    def evaluate(self, present_kos):
        if not self.children:
            return 1.0, True, set(), set()
        scores = []
        all_matched = set()
        all_missing = set()
        all_done = True

        for c in self.children:
            sc, done, matched, missing = c.evaluate(present_kos)
            scores.append(sc)
            all_matched.update(matched)
            all_missing.update(missing)
            if not done:
                all_done = False

        completeness = sum(scores) / len(scores)
        return completeness, all_done, all_matched, all_missing

    def __repr__(self):
        return "(" + " ".join(str(c) for c in self.children) + ")"


def tokenize(definition_str):
    """Tokenize definition string into parentheses, operators, and KO identifiers."""
    # Replace multiple spaces with a single space
    s = re.sub(r'\s+', ' ', definition_str.strip())
    tokens = []
    i = 0
    n = len(s)
    while i < n:
        c = s[i]
        if s[i:i+2] == '--':
            tokens.append('--')
            i += 2
        elif c in ['(', ')', ',', '+', '-']:
            tokens.append(c)
            i += 1
        elif c == ' ':
            tokens.append(' ')
            i += 1
        elif c in ['K', 'k']:
            # KO identifier match: K\d{5}
            m = re.match(r'[Kk]\d{5}', s[i:])
            if m:
                tokens.append(m.group(0).upper())
                i += len(m.group(0))
            else:
                i += 1
        elif c in ['M', 'm']:
            # Module identifier match in Module sets: M\d{5}
            m = re.match(r'[Mm]\d{5}', s[i:])
            if m:
                tokens.append(m.group(0).upper())
                i += len(m.group(0))
            else:
                i += 1
        else:
            i += 1
    return tokens


class ModuleNode(Node):
    def __init__(self, module_id):
        self.module_id = module_id.upper()

    def evaluate(self, present_kos):
        return 0.0, False, set(), {self.module_id}

    def __repr__(self):
        return self.module_id


class EmptyNode(Node):
    def evaluate(self, present_kos):
        return 0.0, False, set(), set()

    def __repr__(self):
        return ""


def parse_definition(definition_str):
    """
    Parses KEGG Module definition grammar with operator precedence:
    Lowest to highest:
    1. Space ' ' (AND across steps)
    2. Comma ',' (OR across alternatives)
    3. Plus '+' / Minus '-' (Complex subunit conjunction / optional subunit)
    4. Minus '-' (Optional prefix)
    5. Parentheses '()'
    """
    tokens = tokenize(definition_str)
    pos = 0

    def parse_and():
        nonlocal pos
        nodes = [parse_or()]
        while pos < len(tokens) and tokens[pos] == ' ':
            pos += 1
            nodes.append(parse_or())
        return nodes[0] if len(nodes) == 1 else AndNode(nodes)

    def parse_or():
        nonlocal pos
        nodes = [parse_complex()]
        while pos < len(tokens) and tokens[pos] == ',':
            pos += 1
            nodes.append(parse_complex())
        return nodes[0] if len(nodes) == 1 else OrNode(nodes)

    def parse_complex():
        nonlocal pos
        first = parse_optional()
        if pos < len(tokens) and tokens[pos] in ('+', '-'):
            children = [first]
            while pos < len(tokens) and tokens[pos] in ('+', '-'):
                op = tokens[pos]
                pos += 1
                if op == '+':
                    children.append(parse_optional())
                elif op == '-':
                    child = parse_atom()
                    children.append(OptionalNode(child))
            return ComplexNode(children)
        return first

    def parse_optional():
        nonlocal pos
        if pos < len(tokens) and tokens[pos] == '-':
            pos += 1
            child = parse_atom()
            return OptionalNode(child)
        return parse_atom()

    def parse_atom():
        nonlocal pos
        if pos >= len(tokens):
            return EmptyNode()
        t = tokens[pos]
        if t == '(':
            pos += 1
            node = parse_and()
            if pos < len(tokens) and tokens[pos] == ')':
                pos += 1
            return node
        elif re.match(r'^K\d{5}$', t):
            pos += 1
            return KONode(t)
        elif re.match(r'^M\d{5}$', t):
            pos += 1
            return ModuleNode(t)
        elif t == '--':
            pos += 1
            return OptionalNode(EmptyNode())
        else:
            pos += 1
            return EmptyNode()

    res = parse_and()
    if pos < len(tokens):
        extra_nodes = [res]
        while pos < len(tokens):
            if tokens[pos] == ' ':
                pos += 1
                continue
            extra_nodes.append(parse_or())
        if len(extra_nodes) > 1:
            res = AndNode(extra_nodes)
    return res


def evaluate_module(definition_str_or_tree, present_kos, partial_score_threshold=0.50, max_missing_steps=1):
    """
    Evaluates a KEGG module definition (or pre-parsed AST Node tree) against a set of present KOs.
    Scientific completeness criteria:
    - complete: definition_satisfied == True (100% of required Boolean blocks satisfied)
    - near_complete: definition_satisfied == False, but missing_required_blocks <= max_missing_steps (e.g. missing <=1 step in >=3-step pathways)
    - partial: completeness_score >= partial_score_threshold (default: 0.50)
    - absent: completeness_score < partial_score_threshold and not near_complete

    Returns:
      completeness (float 0.0-1.0),
      status ("complete", "near_complete", "partial", "absent"),
      definition_satisfied (bool),
      matched_kos (sorted list),
      missing_kos (sorted list)
    """
    if isinstance(definition_str_or_tree, Node):
        tree = definition_str_or_tree
    else:
        tree = parse_definition(definition_str_or_tree)
    completeness, definition_satisfied, matched, missing = tree.evaluate(set(present_kos))
    completeness = round(min(1.0, max(0.0, completeness)), 4)

    missing_one_step = False
    if not definition_satisfied and isinstance(tree, AndNode) and len(tree.children) >= 3 and max_missing_steps > 0:
        completed_steps = sum(1 for c in tree.children if c.evaluate(set(present_kos))[1])
        if len(tree.children) - completed_steps <= max_missing_steps:
            missing_one_step = True

    if definition_satisfied:
        status = "complete"
    elif missing_one_step:
        status = "near_complete"
    elif completeness >= partial_score_threshold:
        status = "partial"
    else:
        status = "absent"

    return completeness, status, definition_satisfied, sorted(matched), sorted(missing)
