
import os
import Bio
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio import Restriction
from Bio.Restriction import *
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqFeature
from Bio.SeqFeature import *
import eating
import sqlite3

def db_format(dbname):
    '''
    dbname: the filename of the sqlite db to store guide scores in. Can be an existing database, in which case
        new entries are appended. If the filename isn't found, a new db with this name is created.
    '''
    # Generate a DictCursor (i.e.
    def dict_factory(cursor, row):
        d = {}
        for idx, col in enumerate(cursor.description):
            d[col[0]] = row[idx]
        return d

    class DictConnnection(sqlite3.Connection):
        def __init__(self, *args, **kwargs):
            sqlite3.Connection.__init__(self, *args, **kwargs)

        def cursor(self):
            return DictCursor(self)

    class DictCursor(sqlite3.Cursor):
        def __init__(self, *args, **kwargs):
            sqlite3.Cursor.__init__(self, *args, **kwargs)
            self.row_factory = lambda cur, row: dict_factory(self, row)

    if os.path.isfile(str(dbname)):
        print("Database " + dbname +" exists. Appending to this database.")
    #
    con = sqlite3.connect(str(dbname), factory=DictConnnection)

    try:
        with con:
            con.row_factory = dict_factory
            cur = con.cursor()
            print("Creating tables.")
            cur.execute('''CREATE TABLE scores
                     (sequence text PRIMARY KEY, genome text, score real)''')
            cur.execute('''CREATE TABLE locations
                     (id integer, sequence text, genome text, loc_start integer, scaffold test, enzyme text, rel_name text)''')
    except:
        print("Tables already set up.")
    return None

def get_score(guide, database):
    con = sqlite3.connect(str(database))
    with con:
        con.row_factory = sqlite3.Row
        c = con.cursor()
        wasnt_list_at_start = 0
        # Operate on lists if that's what's given
        if isinstance(guide, list) == False:
            guide = [guide]
            wasnt_list_at_start = 1
        failed = []
        for item in enumerate(guide):
            try:
                c.execute('SELECT * FROM scores WHERE sequence="{}"'.format(item[1].seq))
            except: # This isn't used yet...
                failed.append(item[0])
        all_rows = c.fetchall()

    output = {}
    # Get the cursor data into a list of tuples
    for row in all_rows:
        output[str(row["sequence"])] = row["score"]

    # Apply the data in the tuples back to the Seq objects originally passed
    for item in guide:
        item.annotations["score"] = output[str(item.seq)]

    # Give back a single guide if that's what we started with
    if wasnt_list_at_start == 1:
        guide = guide[0]

    return(guide)

def put_score(guides, database, genome):
    con = sqlite3.connect(database, check_same_thread=False)
    with con:
        cursor = con.cursor()
        for guide in guides:
            try:
                cursor.execute("INSERT INTO scores(sequence, genome, score) VALUES('{}', '{}', '{}')".format(guide.seq, genome, guide.annotations['score']))
            except sqlite3.Error as e:
                print "An error occurred:", e.args[0], guide.seq

def presence_check(item, sqlitedb):
    present = 0
    try:
        success = eating.sqlite.get_score(item, sqlitedb)
        present = 1
    except:
        present = 0
    return present
