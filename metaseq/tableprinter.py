import fisher
import numpy as np

def print_2x2_table(table, row_labels, col_labels, fmt="%d"):
    """
    Prints a table used for Fisher's exact test. Adds row, column, and grand
    totals.

    :param table: The four cells of a 2x2 table: [r1c1, r1c2, r2c1, r2c2]
    :param row_labels: A length-2 list of row names
    :param col_labels: A length-2 list of column names

    """
    grand = sum(table)

    # Separate table into components and get row/col sums
    t11, t12, t21, t22 = table

    # Row sums, col sums, and grand total
    r1 = t11 + t12
    r2 = t21 + t22
    c1 = t11 + t21
    c2 = t12 + t22

    # Re-cast everything as the appropriate format
    t11, t12, t21, t22, c1, c2, r1, r2, grand = [
        fmt % i for i in [t11, t12, t21, t22, c1, c2, r1, r2, grand]]

    # Construct rows and columns the long way...
    rows = [
        [""] + col_labels + ['total'],
        [row_labels[0], t11, t12, r1],
        [row_labels[1], t21, t22, r2],
        ['total', c1, c2, grand],
    ]

    cols = [
        [row[0] for row in rows],
        [col_labels[0], t11, t21, c1],
        [col_labels[1], t12, t22, c2],
        ['total', r1, r2, grand],
    ]

    # Get max column width for each column; need this for nice justification
    widths = []
    for col in cols:
        widths.append(max(len(i) for i in col))

    # ReST-formatted header
    sep = ['=' * i for i in widths]

    # Construct the table one row at a time with nice justification
    s = []
    s.append(' '.join(sep))
    s.append(' '.join(i.ljust(j) for i, j in zip(rows[0], widths)))
    s.append(' '.join(sep))
    for row in rows[1:]:
        s.append(' '.join(i.ljust(j) for i, j in zip(row, widths)))
    s.append(' '.join(sep) + '\n')
    return "\n".join(s)


def print_row_perc_table(table, row_labels, col_labels):
    """
    given a table, print the percentages rather than the totals
    """
    r1c1, r1c2, r2c1, r2c2 = map(float, table)
    row1 = r1c1 + r1c2
    row2 = r2c1 + r2c2

    blocks = [
        (r1c1, row1),
        (r1c2, row1),
        (r2c1, row2),
        (r2c2, row2)]

    new_table = []

    for cell, row in blocks:
        try:
            x = cell / row
        except ZeroDivisionError:
            x = 0
        new_table.append(x)

    s = print_2x2_table(new_table, row_labels, col_labels, fmt="%.2f")
    s = s.splitlines(True)
    del s[5]
    return ''.join(s)


def print_col_perc_table(table, row_labels, col_labels):
    """
    given a table, print the cols as percentages
    """
    r1c1, r1c2, r2c1, r2c2 = map(float, table)
    col1 = r1c1 + r2c1
    col2 = r1c2 + r2c2


    blocks = [
        (r1c1, col1),
        (r1c2, col2),
        (r2c1, col1),
        (r2c2, col2)]

    new_table = []

    for cell, row in blocks:
        try:
            x = cell / row
        except ZeroDivisionError:
            x = 0
        new_table.append(x)

    s = print_2x2_table(new_table, row_labels, col_labels, fmt="%.2f")
    s = s.splitlines(False)
    last_space = s[0].rindex(" ")
    new_s = [i[:last_space] for i in s]
    return '\n'.join(new_s)


def table_maker(subset, ind1, ind2, row_labels, col_labels, title):
    """
    `subset` provides a subsetted boolean of items to consider.  If no subset,
    you can use all with `np.ones_like(ind1) == 1`

    `ind1` is used to subset rows, e.g., log2fc > 0.  This is used for rows, so
    row_label might be ['upregulated', 'others']

    `ind2` is used to subset cols.  For example, col_labels would be
    ['bound', 'unbound']
    """
    table = [
        sum(subset & ind1 & ind2),
        sum(subset & ind1 & ~ind2),
        sum(subset & ~ind1 & ind2),
        sum(subset & ~ind1 & ~ind2)
    ]
    print
    print title
    print '-' * len(title)
    print print_2x2_table(table, row_labels=row_labels, col_labels=col_labels)
    print print_row_perc_table(
        table, row_labels=row_labels, col_labels=col_labels)
    print print_col_perc_table(
        table, row_labels=row_labels, col_labels=col_labels)
    print fisher.pvalue(*table)


if __name__ == "__main__":
    table = [12, 5, 29, 2]
    s = print_2x2_table(
        table,
        row_labels=['Selected', 'Not selected'],
        col_labels=['Having the property', 'Not having the property']
    )

    str_table = """
    ============ =================== ======================= =====
                 Having the property Not having the property total
    ============ =================== ======================= =====
    Selected     12                  5                       17
    Not selected 29                  2                       31
    total        41                  7                       48
    ============ =================== ======================= =====
    """

    # For the test, remove the first newline and all common leading whitespace
    from textwrap import dedent
    str_table = "".join(str_table.splitlines(True)[1:])
    print s
    assert dedent(str_table) == s
