#!/usr/bin python3
"""
#
#------------------------------------------------------------------------------
#                         University of Minnesota
#         Copyright 2017, Regents of the University of Minnesota
#------------------------------------------------------------------------------
# Author:
#
#  James E Johnson
#
#------------------------------------------------------------------------------
"""

from __future__ import print_function
from __future__ import unicode_literals

import argparse
import re
import sys

from Bio.Seq import translate
from Bio.Seq import reverse_complement

import digest

from twobitreader import TwoBitFile

from time import sleep
import requests


server = "https://rest.ensembl.org"
ext = "/info/assembly/homo_sapiens?"
max_region = 4000000
debug = False


import itertools as it
from collections import deque


def cleave(sequence, rule, missed_cleavages=0, min_length=None):
    """Cleaves a polypeptide sequence using a given rule.

    Parameters
    ----------
    sequence : str
        The sequence of a polypeptide.

        .. note::
            The sequence is expected to be in one-letter uppercase notation.
            Otherwise, some of the cleavage rules in :py:data:`expasy_rules`
            will not work as expected.

    rule : str or compiled regex
        A regular expression describing the site of cleavage. It is recommended
        to design the regex so that it matches only the residue whose
        C-terminal bond is to be cleaved. All additional requirements should be
        specified using `lookaround assertions
        <http://www.regular-expressions.info/lookaround.html>`_.
        :py:data:`expasy_rules` contains cleavage rules
        for popular cleavage agents.
    missed_cleavages : int, optional
        Maximum number of allowed missed cleavages. Defaults to 0.
    min_length : int or None, optional
        Minimum peptide length. Defaults to :py:const:`None`.

        ..note ::
            This checks for string length, which is only correct for one-letter
            notation and not for full *modX*. Use :py:func:`length` manually if
            you know what you are doing and apply :py:func:`cleave` to *modX*
            sequences.

    Returns
    -------
    out : set
        A set of unique (!) peptides.

    Examples
    --------
    >>> cleave('AKAKBK', expasy_rules['trypsin'], 0) == {'AK', 'BK'}
    True
    >>> cleave('GKGKYKCK', expasy_rules['trypsin'], 2) == \
    {'CK', 'GKYK', 'YKCK', 'GKGK', 'GKYKCK', 'GK', 'GKGKYK', 'YK'}
    True

    """
    return set(_cleave(sequence, rule, missed_cleavages, min_length))


def _cleave(sequence, rule, missed_cleavages=0, min_length=None):
    """Like :py:func:`cleave`, but the result is a list. Refer to
    :py:func:`cleave` for explanation of parameters.
    """
    peptides = []
    ml = missed_cleavages+2
    trange = range(ml)
    cleavage_sites = deque([0], maxlen=ml)
    cl = 1
    for i in it.chain([x.end() for x in re.finditer(rule, sequence)],
                      [None]):
        cleavage_sites.append(i)
        if cl < ml:
            cl += 1
        for j in trange[:cl-1]:
            seq = sequence[cleavage_sites[j]:cleavage_sites[-1]]
            if seq:
                if min_length is None or len(seq) >= min_length:
                    peptides.append(seq)
    return peptides


def num_sites(sequence, rule, **kwargs):
    """Count the number of sites where `sequence` can be cleaved using
    the given `rule` (e.g. number of miscleavages for a peptide).

    Parameters
    ----------
    sequence : str
        The sequence of a polypeptide.
    rule : str or compiled regex
        A regular expression describing the site of cleavage. It is recommended
        to design the regex so that it matches only the residue whose
        C-terminal bond is to be cleaved. All additional requirements should be
        specified using `lookaround assertions
        <http://www.regular-expressions.info/lookaround.html>`_.
    labels : list, optional
        A list of allowed labels for amino acids and terminal modifications.

    Returns
    -------
    out : int
        Number of cleavage sites.
    """
    return len(_cleave(sequence, rule, **kwargs)) - 1


expasy_rules = {
    'arg-c':         r'R',
    'asp-n':         r'\w(?=D)',
    'bnps-skatole': r'W',
    'caspase 1':     r'(?<=[FWYL]\w[HAT])D(?=[^PEDQKR])',
    'caspase 2':     r'(?<=DVA)D(?=[^PEDQKR])',
    'caspase 3':     r'(?<=DMQ)D(?=[^PEDQKR])',
    'caspase 4':     r'(?<=LEV)D(?=[^PEDQKR])',
    'caspase 5':     r'(?<=[LW]EH)D',
    'caspase 6':     r'(?<=VE[HI])D(?=[^PEDQKR])',
    'caspase 7':     r'(?<=DEV)D(?=[^PEDQKR])',
    'caspase 8':     r'(?<=[IL]ET)D(?=[^PEDQKR])',
    'caspase 9':     r'(?<=LEH)D',
    'caspase 10':    r'(?<=IEA)D',
    'chymotrypsin high specificity': r'([FY](?=[^P]))|(W(?=[^MP]))',
    'chymotrypsin low specificity':
        r'([FLY](?=[^P]))|(W(?=[^MP]))|(M(?=[^PY]))|(H(?=[^DMPW]))',
    'clostripain':   r'R',
    'cnbr':          r'M',
    'enterokinase':  r'(?<=[DE]{3})K',
    'factor xa':     r'(?<=[AFGILTVM][DE]G)R',
    'formic acid':   r'D',
    'glutamyl endopeptidase': r'E',
    'granzyme b':    r'(?<=IEP)D',
    'hydroxylamine': r'N(?=G)',
    'iodosobenzoic acid': r'W',
    'lysc':          r'K',
    'ntcb':          r'\w(?=C)',
    'pepsin ph1.3':  r'((?<=[^HKR][^P])[^R](?=[FLWY][^P]))|'
                     r'((?<=[^HKR][^P])[FLWY](?=\w[^P]))',
    'pepsin ph2.0':  r'((?<=[^HKR][^P])[^R](?=[FL][^P]))|'
                     r'((?<=[^HKR][^P])[FL](?=\w[^P]))',
    'proline endopeptidase': r'(?<=[HKR])P(?=[^P])',
    'proteinase k':  r'[AEFILTVWY]',
    'staphylococcal peptidase i': r'(?<=[^E])E',
    'thermolysin':   r'[^DE](?=[AFILMV])',
    'thrombin':      r'((?<=G)R(?=G))|'
                     r'((?<=[AFGILTVM][AFGILTVWA]P)R(?=[^DE][^DE]))',
    'trypsin':       r'([KR](?=[^P]))|((?<=W)K(?=P))|((?<=M)R(?=P))'
    }
"""
This dict contains regular expressions for cleavage rules of the most
popular proteolytic enzymes. The rules were taken from the
`PeptideCutter tool
<http://ca.expasy.org/tools/peptidecutter/peptidecutter_enzymes.html>`_
at Expasy.
"""


def ensembl_rest(ext, headers):
    if debug:
        print("%s" % ext, file=sys.stderr)
    r = requests.get(server+ext, headers=headers)
    if r.status_code == 429:
        print("response headers: %s\n" % r.headers, file=sys.stderr)
        if 'Retry-After' in r.headers:
            sleep(r.headers['Retry-After'])
            r = requests.get(server+ext, headers=headers)
    if not r.ok:
        r.raise_for_status()
    return r


def get_species():
    results = dict()
    ext = "/info/species"
    req_header = {"Content-Type": "application/json"}
    r = ensembl_rest(ext, req_header)
    for species in r.json()['species']:
        results[species['name']] = species
        print("%s\t%s\t%s\t%s\t%s" %
              (species['name'], species['common_name'],
               species['display_name'],
               species['strain'],
               species['taxon_id']), file=sys.stdout)
    return results


def get_biotypes(species):
    biotypes = []
    ext = "/info/biotypes/%s?" % species
    req_header = {"Content-Type": "application/json"}
    r = ensembl_rest(ext, req_header)
    for entry in r.json():
        if 'biotype' in entry:
            biotypes.append(entry['biotype'])
    return biotypes


def get_toplevel(species):
    coord_systems = dict()
    ext = "/info/assembly/%s?" % species
    req_header = {"Content-Type": "application/json"}
    r = ensembl_rest(ext, req_header)
    toplevel = r.json()
    for seq in toplevel['top_level_region']:
        if seq['coord_system'] not in coord_systems:
            coord_systems[seq['coord_system']] = dict()
        coord_system = coord_systems[seq['coord_system']]
        coord_system[seq['name']] = int(seq['length'])
    return coord_systems


def get_transcripts_bed(species, refseq, start, length, strand='',
                        params=None):
    bed = []
    param = params if params else ''
    req_header = {"Content-Type": "text/x-bed"}
    regions = list(range(start, length, max_region))
    if not regions or regions[-1] < length:
        regions.append(length)
    for end in regions[1:]:
        ext = "/overlap/region/%s/%s:%d-%d%s?feature=transcript;%s"\
            % (species, refseq, start, end, strand, param)
        start = end + 1
        r = ensembl_rest(ext, req_header)
        if r.text:
            bed += r.text.splitlines()
    return bed


def get_seq(id, seqtype, params=None):
    param = params if params else ''
    ext = "/sequence/id/%s?type=%s;%s" % (id, seqtype, param)
    req_header = {"Content-Type": "text/plain"}
    r = ensembl_rest(ext, req_header)
    return r.text


def get_cdna(id, params=None):
    return get_seq(id, 'cdna', params=params)


def get_cds(id, params=None):
    return get_seq(id, 'cds', params=params)


def get_genomic(id, params=None):
    return get_seq(id, 'genomic', params=params)


def get_transcript_haplotypes(species, transcript):
    ext = "/transcript_haplotypes/%s/%s?aligned_sequences=1"\
        % (species, transcript)
    req_header = {"Content-Type": "application/json"}
    r = ensembl_rest(ext, req_header)
    decoded = r.json()
    return decoded

def bed_from_line(line, ensembl=False, seq_column=None):
    fields = line.rstrip('\r\n').split('\t')
    if len(fields) < 12:
        return None
    (chrom, chromStart, chromEnd, name, score, strand,
     thickStart, thickEnd, itemRgb,
     blockCount, blockSizes, blockStarts) = fields[0:12]
    bed_entry = BedEntry(chrom=chrom, chromStart=chromStart, chromEnd=chromEnd,
                         name=name, score=score, strand=strand,
                         thickStart=thickStart, thickEnd=thickEnd,
                         itemRgb=itemRgb,
                         blockCount=blockCount,
                         blockSizes=blockSizes.rstrip(','),
                         blockStarts=blockStarts.rstrip(','))
    if seq_column is not None and -len(fields) <= seq_column < len(fields):
        bed_entry.seq = fields[seq_column]
    if ensembl and len(fields) >= 20:
        bed_entry.second_name = fields[12]
        bed_entry.cds_start_status = fields[13]
        bed_entry.cds_end_status = fields[14]
        bed_entry.exon_frames = fields[15].rstrip(',')
        bed_entry.biotype = fields[16]
        bed_entry.gene_name = fields[17]
        bed_entry.second_gene_name = fields[18]
        bed_entry.gene_type = fields[19]
    return bed_entry


def as_int_list(obj):
    if obj is None:
        return None
    if isinstance(obj, list):
        return [int(x) for x in obj]
    elif isinstance(obj, str):
        return [int(x) for x in obj.split(',')]
    else:  # python2 unicode?
        return [int(x) for x in str(obj).split(',')]


class BedEntry(object):
    def __init__(self, chrom=None, chromStart=None, chromEnd=None,
                 name=None, score=None, strand=None,
                 thickStart=None, thickEnd=None, itemRgb=None,
                 blockCount=None, blockSizes=None, blockStarts=None):
        self.chrom = chrom
        self.chromStart = int(chromStart)
        self.chromEnd = int(chromEnd)
        self.name = name
        self.score = int(score) if score is not None else 0
        self.strand = '-' if str(strand).startswith('-') else '+'
        self.thickStart = int(thickStart) if thickStart else self.chromStart
        self.thickEnd = int(thickEnd) if thickEnd else self.chromEnd
        self.itemRgb = str(itemRgb) if itemRgb is not None else r'100,100,100'
        self.blockCount = int(blockCount)
        self.blockSizes = as_int_list(blockSizes)
        self.blockStarts = as_int_list(blockStarts)
        self.second_name = None
        self.cds_start_status = None
        self.cds_end_status = None
        self.exon_frames = None
        self.biotype = None
        self.gene_name = None
        self.second_gene_name = None
        self.gene_type = None
        self.seq = None
        self.cdna = None
        self.pep = None
        # T26C
        self.aa_change = []
        # p.Trp26Cys g.<pos><ref>><alt> # g.1304573A>G
        self.variants = []

    def __str__(self):
        return '%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%s\t%d\t%s\t%s' % (
            self.chrom, self.chromStart, self.chromEnd,
            self.name, self.score, self.strand,
            self.thickStart, self.thickEnd, str(self.itemRgb), self.blockCount,
            ','.join([str(x) for x in self.blockSizes]),
            ','.join([str(x) for x in self.blockStarts]))

    def get_splice_junctions(self):
        splice_juncs = []
        for i in range(self.blockCount - 1):
            splice_junc = "%s:%d_%d"\
                % (self.chrom,
                   self.chromStart + self.blockSizes[i],
                   self.chromStart + self.blockStarts[i+1])
            splice_juncs.append(splice_junc)
        return splice_juncs

    def get_exon_seqs(self):
        if not self.seq:
            return None
        exons = []
        for i in range(self.blockCount):
            exons.append(self.seq[self.blockStarts[i]:self.blockStarts[i]
                         + self.blockSizes[i]])
        if self.strand == '-':  # reverse complement
            exons.reverse()
            for i, s in enumerate(exons):
                exons[i] = reverse_complement(s)
        return exons

    def get_spliced_seq(self, strand=None):
        if not self.seq:
            return None
        seq = ''.join(self.get_exon_seqs())
        if strand and self.strand != strand:
            seq = reverse_complement(seq)
        return seq

    def get_cdna(self):
        if not self.cdna:
            self.cdna = self.get_spliced_seq()
        return self.cdna

    def get_cds(self):
        cdna = self.get_cdna()
        if cdna:
            if self.chromStart == self.thickStart\
               and self.chromEnd == self.thickEnd:
                return cdna
            pos = [self.cdna_offset_of_pos(self.thickStart),
                   self.cdna_offset_of_pos(self.thickEnd)]
            if 0 <= min(pos) <= max(pos) <= len(cdna):
                return cdna[min(pos):max(pos)]
        return None

    def set_cds(self, cdna_start, cdna_end):
        cdna_len = sum(self.blockSizes)
        if 0 <= cdna_start < cdna_end <= cdna_len:
            cds_pos = [self.pos_of_cdna_offet(cdna_start),
                       self.pos_of_cdna_offet(cdna_end)]
            if all(cds_pos):
                self.thickStart = min(cds_pos)
                self.thickEnd = max(cds_pos)
            return self
        return None

    def trim_cds(self, basepairs):
        if self.chromStart <= self.thickStart < self.thickEnd <= self.chromEnd:
            cds_pos = [self.cdna_offset_of_pos(self.thickStart),
                       self.cdna_offset_of_pos(self.thickEnd)]
            if basepairs > 0:
                return self.set_cds(min(cds_pos) + basepairs, max(cds_pos))
            else:
                return self.set_cds(min(cds_pos), max(cds_pos) + basepairs)
        return None

    def get_cds_bed(self):
        cds_pos = [self.cdna_offset_of_pos(self.thickStart),
                   self.cdna_offset_of_pos(self.thickEnd)]
        return self.trim(min(cds_pos), max(cds_pos))

    def get_cigar(self):
        cigar = ''
        r = range(self.blockCount)
        xl = None
        for x in r:
            if xl is not None:
                intronSize = abs(self.blockStarts[x] - self.blockSizes[xl]
                                 - self.blockStarts[xl])
                cigar += '%dN' % intronSize
            cigar += '%dM' % self.blockSizes[x]
            xl = x
        return cigar

    def get_cigar_md(self):
        cigar = ''
        md = ''
        r = range(self.blockCount)
        xl = None
        for x in r:
            if xl is not None:
                intronSize = abs(self.blockStarts[x] - self.blockSizes[xl]
                                 - self.blockStarts[xl])
                cigar += '%dN' % intronSize
            cigar += '%dM' % self.blockSizes[x]
            xl = x
        md = '%d' % sum(self.blockSizes)
        return (cigar, md)

    def get_translation(self, sequence=None):
        translation = None
        seq = sequence if sequence else self.get_spliced_seq()
        if seq:
            seqlen = len(seq) / 3 * 3
            if seqlen >= 3:
                translation = translate(seq[:seqlen])
        return translation

    def get_translations(self):
        translations = []
        seq = self.get_spliced_seq()
        if seq:
            for i in range(3):
                translation = self.get_translation(sequence=seq[i:])
                if translation:
                    translations.append(translation)
        return translations

    def pos_of_cdna_offet(self, offset):
        if offset is not None and 0 <= offset < sum(self.blockSizes):
            r = list(range(self.blockCount))
            rev = self.strand == '-'
            if rev:
                r.reverse()
            nlen = 0
            for x in r:
                if offset < nlen + self.blockSizes[x]:
                    if rev:
                        return self.chromStart + self.blockStarts[x]\
                               + self.blockSizes[x] - (offset - nlen)
                    else:
                        return self.chromStart + self.blockStarts[x]\
                               + (offset - nlen)
                nlen += self.blockSizes[x]
        return None

    def cdna_offset_of_pos(self, pos):
        if not self.chromStart <= pos < self.chromEnd:
            return -1
        r = list(range(self.blockCount))
        rev = self.strand == '-'
        if rev:
            r.reverse()
        nlen = 0
        for x in r:
            bStart = self.chromStart + self.blockStarts[x]
            bEnd = bStart + self.blockSizes[x]
            if bStart <= pos < bEnd:
                return nlen + (bEnd - pos if rev else pos - bStart)
            nlen += self.blockSizes[x]

    def apply_variant(self, pos, ref, alt):
        pos = int(pos)
        if not ref or not alt:
            print("variant requires ref and alt sequences", file=sys.stderr)
            return
        if not self.chromStart <= pos <= self.chromEnd:
            print("variant not in entry %s: %s %d < %d < %d" %
                  (self.name, self.strand,
                   self.chromStart, pos, self.chromEnd),
                  file=sys.stderr)
            print("%s" % str(self), file=sys.stderr)
            return
        if len(ref) != len(alt):
            print("variant only works for snp: %s  %s" % (ref, alt),
                  file=sys.stderr)
            return
        if not self.seq:
            print("variant entry %s has no seq" % self.name, file=sys.stderr)
            return
        """
        if self.strand  == '-':
            ref = reverse_complement(ref)
            alt = reverse_complement(alt)
        """
        bases = list(self.seq)
        offset = pos - self.chromStart
        for i in range(len(ref)):
            # offset = self.cdna_offset_of_pos(pos+i)
            if offset is not None:
                bases[offset+i] = alt[i]
            else:
                print("variant offset %s: %s %d < %d < %d" %
                      (self.name, self.strand, self.chromStart,
                       pos+1, self.chromEnd), file=sys.stderr)
                print("%s" % str(self), file=sys.stderr)
        self.seq = ''.join(bases)
        self.variants.append("g.%d%s>%s" % (pos+1, ref, alt))

    def get_variant_bed(self, pos, ref, alt):
        pos = int(pos)
        if not ref or not alt:
            print("variant requires ref and alt sequences", file=sys.stderr)
            return None
        if not self.chromStart <= pos <= self.chromEnd:
            print("variant not in entry %s: %s %d < %d < %d" %
                  (self.name, self.strand,
                   self.chromStart, pos, self.chromEnd),
                  file=sys.stderr)
            print("%s" % str(self), file=sys.stderr)
            return None
        if not self.seq:
            print("variant entry %s has no seq" % self.name, file=sys.stderr)
            return None
        tbed = BedEntry(chrom=self.chrom,
                        chromStart=self.chromStart, chromEnd=self.chromEnd,
                        name=self.name, score=self.score, strand=self.strand,
                        thickStart=self.chromStart, thickEnd=self.chromEnd,
                        itemRgb=self.itemRgb,
                        blockCount=self.blockCount,
                        blockSizes=self.blockSizes,
                        blockStarts=self.blockStarts)
        bases = list(self.seq)
        offset = pos - self.chromStart
        tbed.seq = ''.join(bases[:offset] + list(alt)
                           + bases[offset+len(ref):])
        if len(ref) != len(alt):
            diff = len(alt) - len(ref)
            rEnd = pos + len(ref)
            # need to adjust blocks
            #  change spans blocks,
            for x in range(tbed.blockCount):
                bStart = tbed.chromStart + tbed.blockStarts[x]
                bEnd = bStart + tbed.blockSizes[x]
                # change within a block or extends (last block)
                #  adjust blocksize
                #  seq:            GGGcatGGG
                #  ref c alt tag:  GGGtagatGGG
                #  ref cat alt a:  GGGaGGG
                if bStart <= pos < rEnd < bEnd:
                    tbed.blockSizes[x] += diff
        return tbed

    # (start, end)
    def get_subrange(self, tstart, tstop, debug=False):
        chromStart = self.chromStart
        chromEnd = self.chromEnd
        if debug:
            print("%s" % (str(self)), file=sys.stderr)
        r = list(range(self.blockCount))
        if self.strand == '-':
            r.reverse()
        bStart = 0
        bEnd = 0
        for x in r:
            bEnd = bStart + self.blockSizes[x]
            if bStart <= tstart < bEnd:
                if self.strand == '+':
                    chromStart = self.chromStart + self.blockStarts[x] +\
                        (tstart - bStart)
                else:
                    chromEnd = self.chromStart + self.blockStarts[x] +\
                        self.blockSizes[x] - (tstart - bStart)
            if bStart <= tstop < bEnd:
                if self.strand == '+':
                    chromEnd = self.chromStart + self.blockStarts[x] +\
                        (tstop - bStart)
                else:
                    chromStart = self.chromStart + self.blockStarts[x] +\
                        self.blockSizes[x] - (tstop - bStart)
            if debug:
                print("%3d %s\t%d\t%d\t%d\t%d\t%d\t%d" %
                      (x, self.strand, bStart, bEnd,
                       tstart, tstop, chromStart, chromEnd), file=sys.stderr)
            bStart += self.blockSizes[x]
        return(chromStart, chromEnd)

    # get the blocks for sub range
    def get_blocks(self, chromStart, chromEnd):
        tblockCount = 0
        tblockSizes = []
        tblockStarts = []
        for x in range(self.blockCount):
            bStart = self.chromStart + self.blockStarts[x]
            bEnd = bStart + self.blockSizes[x]
            if bStart > chromEnd:
                break
            if bEnd < chromStart:
                continue
            cStart = max(chromStart, bStart)
            tblockStarts.append(cStart - chromStart)
            tblockSizes.append(min(chromEnd, bEnd) - cStart)
            tblockCount += 1
        return (tblockCount, tblockSizes, tblockStarts)

    def trim(self, tstart, tstop, debug=False):
        (tchromStart, tchromEnd) =\
            self.get_subrange(tstart, tstop, debug=debug)
        (tblockCount, tblockSizes, tblockStarts) =\
            self.get_blocks(tchromStart, tchromEnd)
        tbed = BedEntry(
            chrom=self.chrom, chromStart=tchromStart, chromEnd=tchromEnd,
            name=self.name, score=self.score, strand=self.strand,
            thickStart=tchromStart, thickEnd=tchromEnd, itemRgb=self.itemRgb,
            blockCount=tblockCount,
            blockSizes=tblockSizes, blockStarts=tblockStarts)
        if self.seq:
            ts = tchromStart-self.chromStart
            te = tchromEnd - tchromStart + ts
            tbed.seq = self.seq[ts:te]
        return tbed

    def get_filtered_translations(self, untrimmed=False, filtering=True,
                                  ignore_left_bp=0, ignore_right_bp=0,
                                  debug=False):
        translations = [None, None, None]
        seq = self.get_spliced_seq()
        ignore = (ignore_left_bp if self.strand == '+'
                  else ignore_right_bp) / 3
        block_sum = sum(self.blockSizes)
        exon_sizes = [x for x in self.blockSizes]
        if self.strand == '-':
            exon_sizes.reverse()
        splice_sites = [sum(exon_sizes[:x]) / 3
                        for x in range(1, len(exon_sizes))]
        if debug:
            print("splice_sites: %s" % splice_sites, file=sys.stderr)
        junc = splice_sites[0] if len(splice_sites) > 0 else exon_sizes[0]
        if seq:
            for i in range(3):
                translation = self.get_translation(sequence=seq[i:])
                if translation:
                    tstart = 0
                    tstop = len(translation)
                    offset = (block_sum - i) % 3
                    if debug:
                        print("frame: %d\ttstart: %d  tstop: %d  " +
                              "offset: %d\t%s" %
                              (i, tstart, tstop, offset, translation),
                              file=sys.stderr)
                    if not untrimmed:
                        tstart = translation.rfind('*', 0, junc) + 1
                        stop = translation.find('*', junc)
                        tstop = stop if stop >= 0 else len(translation)
                    offset = (block_sum - i) % 3
                    trimmed = translation[tstart:tstop]
                    if debug:
                        print("frame: %d\ttstart: %d  tstop: %d  " +
                              "offset: %d\t%s" %
                              (i, tstart, tstop, offset, trimmed),
                              file=sys.stderr)
                    if filtering and tstart > ignore:
                        continue
                    # get genomic locations for start and end
                    if self.strand == '+':
                        chromStart = self.chromStart + i + (tstart * 3)
                        chromEnd = self.chromEnd - offset\
                            - (len(translation) - tstop) * 3
                    else:
                        chromStart = self.chromStart + offset\
                            + (len(translation) - tstop) * 3
                        chromEnd = self.chromEnd - i - (tstart * 3)
                    # get the blocks for this translation
                    (tblockCount, tblockSizes, tblockStarts) =\
                        self.get_blocks(chromStart, chromEnd)
                    translations[i] = (chromStart, chromEnd, trimmed,
                                       tblockCount, tblockSizes, tblockStarts)
                    if debug:
                        print("tblockCount: %d tblockStarts: %s " +
                              "tblockSizes: %s" %
                              (tblockCount, tblockStarts, tblockSizes),
                              file=sys.stderr)
        return translations

    def get_seq_id(self, seqtype='unk:unk', reference='', frame=None):
        # Ensembl fasta ID format
        # >ID SEQTYPE:STATUS LOCATION GENE TRANSCRIPT
        # >ENSP00000328693 pep:splice chromosome:NCBI35:1:904515:910768:1\
        #   gene:ENSG00000158815:transcript:ENST00000328693\
        #    gene_biotype:protein_coding transcript_biotype:protein_coding
        frame_name = ''
        chromStart = self.chromStart
        chromEnd = self.chromEnd
        strand = 1 if self.strand == '+' else -1
        if frame is not None:
            block_sum = sum(self.blockSizes)
            offset = (block_sum - frame) % 3
            frame_name = '_' + str(frame + 1)
            if self.strand == '+':
                chromStart += frame
                chromEnd -= offset
            else:
                chromStart += offset
                chromEnd -= frame
        location = "chromosome:%s:%s:%s:%s:%s"\
            % (reference, self.chrom, chromStart, chromEnd, strand)
        seq_id = "%s%s %s %s" % (self.name, frame_name, seqtype, location)
        return seq_id

    def get_line(self, start_offset=0, end_offset=0):
        if start_offset or end_offset:
            s_offset = start_offset if start_offset else 0
            e_offset = end_offset if end_offset else 0
            if s_offset > self.chromStart:
                s_offset = self.chromStart
            chrStart = self.chromStart - s_offset
            chrEnd = self.chromEnd + e_offset
            blkSizes = self.blockSizes
            blkSizes[0] += s_offset
            blkSizes[-1] += e_offset
            blkStarts = self.blockStarts
            for i in range(1, self.blockCount):
                blkStarts[i] += s_offset
            items = [str(x) for x in [self.chrom, chrStart, chrEnd, self.name,
                                      self.score, self.strand, self.thickStart,
                                      self.thickEnd, self.itemRgb,
                                      self.blockCount,
                                      ','.join([str(x) for x in blkSizes]),
                                      ','.join([str(x) for x in blkStarts])]]
            return '\t'.join(items) + '\n'
        return self.line

def __main__():
    parser = argparse.ArgumentParser(
        description='Translate from BED')
    parser.add_argument(
        'input_bed', default=None,
        help="BED to translate,  '-' for stdin")
    pg_seq = parser.add_argument_group('Genomic sequence source')
    pg_seq.add_argument(
        '-t', '--twobit', default=None,
        help='Genome reference sequence in 2bit format')
    pg_seq.add_argument(
        '-c', '--column', type=int, default=None,
        help='Column offset containing genomic sequence' +
             'between start and stop (-1) for last column')
    pg_out = parser.add_argument_group('Output options')
    pg_out.add_argument(
        '-f', '--fasta', default=None,
        help='Path to output translations.fasta')
    pg_out.add_argument(
        '-b', '--bed', default=None,
        help='Path to output translations.bed')
    pg_bed = parser.add_argument_group('BED filter options')
    pg_bed.add_argument(
        '-E', '--ensembl', action='store_true', default=False,
        help='Input BED is in 20 column Ensembl format')
    pg_bed.add_argument(
        '-R', '--regions', action='append', default=[],
        help='Filter input by regions e.g.:'
             + ' X,2:20000-25000,3:100-500+')
    pg_bed.add_argument(
        '-B', '--biotypes', action='append', default=[],
        help='For Ensembl BED restrict translations to Ensembl biotypes')
    pg_trans = parser.add_argument_group('Translation filter options')
    pg_trans.add_argument(
        '-m', '--min_length', type=int, default=10,
        help='Minimum length of protein translation to report')
    pg_trans.add_argument(
        '-e', '--enzyme', default=None,
        help='Digest translation with enzyme')
    pg_trans.add_argument(
        '-M', '--start_codon', action='store_true', default=False,
        help='Trim translations to methionine start_codon')
    pg_trans.add_argument(
        '-C', '--cds', action='store_true', default=False,
        help='Only translate CDS')
    pg_trans.add_argument(
        '-A', '--all', action='store_true',
        help='Include CDS protein translations ')
    pg_fmt = parser.add_argument_group('ID format options')
    pg_fmt.add_argument(
        '-r', '--reference', default='',
        help='Genome Reference Name')
    pg_fmt.add_argument(
        '-D', '--fa_db', dest='fa_db', default=None,
        help='Prefix DB identifier for fasta ID line, e.g. generic')
    pg_fmt.add_argument(
        '-s', '--fa_sep', dest='fa_sep', default='|',
        help='fasta ID separator defaults to pipe char, ' +
             'e.g. generic|ProtID|description')
    pg_fmt.add_argument(
        '-P', '--id_prefix', default='',
        help='prefix for the sequence ID')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose')
    parser.add_argument('-d', '--debug', action='store_true', help='Debug')
    args = parser.parse_args()

    input_rdr = open(args.input_bed, 'r')\
        if args.input_bed != '-' else sys.stdin
    fa_wtr = open(args.fasta, 'w')\
        if args.fasta is not None and args.fasta != '-' else sys.stdout
    bed_wtr = open(args.bed, 'w') if args.bed is not None else None

    enzyme = expasy_rules.get(args.enzyme, None)

    biotypea = [bt.strip() for biotype in args.biotypes
                for bt in biotype.split(',')]

    twobit = TwoBitFile(args.twobit) if args.twobit else None

    selected_regions = dict()  # chrom:(start, end)
    region_pat = '^(?:chr)?([^:]+)(?::(\d*)(?:-(\d+)([+-])?)?)?'
    if args.regions:
        for entry in args.regions:
            if not entry:
                continue
            regs = [x.strip() for x in entry.split(',') if x.strip()]
            for reg in regs:
                m = re.match(region_pat, reg)
                if m:
                    (chrom, start, end, strand) = m.groups()
                    if chrom:
                        if chrom not in selected_regions:
                            selected_regions[chrom] = []
                        selected_regions[chrom].append([start, end, strand])
        if args.debug:
            print("selected_regions: %s" % selected_regions, file=sys.stderr)

    def filter_by_regions(bed):
        if not selected_regions:
            return True
        ref = re.sub('^(?i)chr', '', bed.chrom)
        if ref not in selected_regions:
            return False
        for reg in selected_regions[ref]:
            (_start, _stop, _strand) = reg
            start = int(_start) if _start else 0
            stop = int(_stop) if _stop else sys.maxint
            if _strand and bed.strand != _strand:
                continue
            if bed.chromEnd >= start and bed.chromStart <= stop:
                return True
        return False

    translations = dict()  # start : end : seq

    def unique_prot(tbed, seq):
        if tbed.chromStart not in translations:
            translations[tbed.chromStart] = dict()
            translations[tbed.chromStart][tbed.chromEnd] = []
            translations[tbed.chromStart][tbed.chromEnd].append(seq)
        elif tbed.chromEnd not in translations[tbed.chromStart]:
            translations[tbed.chromStart][tbed.chromEnd] = []
            translations[tbed.chromStart][tbed.chromEnd].append(seq)
        elif seq not in translations[tbed.chromStart][tbed.chromEnd]:
            translations[tbed.chromStart][tbed.chromEnd].append(seq)
        else:
            return False
        return True

    def get_sequence(chrom, start, end):
        if twobit:
            if chrom in twobit and 0 <= start < end < len(twobit[chrom]):
                return twobit[chrom][start:end]
            contig = chrom[3:] if chrom.startswith('chr') else 'chr%s' % chrom
            if contig in twobit and 0 <= start < end < len(twobit[contig]):
                return twobit[contig][start:end]
        return None

    def write_translation(tbed, accession, peptide):
        if args.id_prefix:
            tbed.name = "%s%s" % (args.id_prefix, tbed.name)
        probed = "%s\t%s\t%s\t%s%s" % (accession, peptide,
                                       'unique', args.reference,
                                       '\t.' * 9)
        if bed_wtr:
            bed_wtr.write("%s\t%s\n" % (str(tbed), probed))
            bed_wtr.flush()
        location = "chromosome:%s:%s:%s:%s:%s"\
            % (args.reference, tbed.chrom,
               tbed.thickStart, tbed.thickEnd, tbed.strand)
        fa_desc = '%s%s' % (args.fa_sep, location)
        fa_db = '%s%s' % (args.fa_db, args.fa_sep) if args.fa_db else ''
        fa_id = ">%s%s%s\n" % (fa_db, tbed.name, fa_desc)
        fa_wtr.write(fa_id)
        fa_wtr.write(peptide)
        fa_wtr.write("\n")
        fa_wtr.flush()

    def translate_bed(bed):
        translate_count = 0
        transcript_id = bed.name
        refprot = None
        if not bed.seq:
            if twobit:
                bed.seq = get_sequence(bed.chrom, bed.chromStart, bed.chromEnd)
            else:
                bed.cdna = get_cdna(transcript_id)
        cdna = bed.get_cdna()
        if cdna is not None:
            cdna_len = len(cdna)
            if args.cds or args.all:
                try:
                    cds = bed.get_cds()
                    if cds:
                        if args.debug:
                            print("cdna:%s" % str(cdna), file=sys.stderr)
                            print("cds: %s" % str(cds), file=sys.stderr)
                        if len(cds) % 3 != 0:
                            cds = cds[:-(len(cds) % 3)]
                        refprot = translate(cds) if cds else None
                except:
                    refprot = None
                if args.cds:
                    if refprot:
                        tbed = bed.get_cds_bed()
                        if args.start_codon:
                            m = refprot.find('M')
                            if m < 0:
                                return 0
                            elif m > 0:
                                bed.trim_cds(m*3)
                                refprot = refprot[m:]
                        stop = refprot.find('*')
                        if stop >= 0:
                            bed.trim_cds((stop - len(refprot)) * 3)
                            refprot = refprot[:stop]
                        if len(refprot) >= args.min_length:
                            write_translation(tbed, bed.name, refprot)
                            return 1
                    return 0
            if args.debug:
                print("%s\n" % (str(bed)), file=sys.stderr)
                print("CDS: %s %d %d" %
                      (bed.strand, bed.cdna_offset_of_pos(bed.thickStart),
                       bed.cdna_offset_of_pos(bed.thickEnd)),
                      file=sys.stderr)
                print("refprot: %s" % str(refprot), file=sys.stderr)
            for offset in range(1):
                seqend = cdna_len - (cdna_len - offset) % 3
                aaseq = translate(cdna[offset:seqend])
                aa_start = 0
                while aa_start < len(aaseq):
                    aa_end = aaseq.find('*', aa_start)
                    if aa_end < 0:
                        aa_end = len(aaseq)
                    prot = aaseq[aa_start:aa_end]
                    if args.start_codon:
                        m = prot.find('M')
                        aa_start += m if m >= 0 else aa_end
                        prot = aaseq[aa_start:aa_end]
                    if enzyme and refprot:
                        frags = _cleave(prot, enzyme)
                        for frag in reversed(frags):
                            if frag in refprot:
                                prot = prot[:prot.rfind(frag)]
                            else:
                                break
                    is_cds = refprot and prot in refprot
                    if args.debug:
                        print("is_cds: %s %s" % (str(is_cds), str(prot)),
                              file=sys.stderr)
                    if len(prot) < args.min_length:
                        pass
                    elif not args.all and is_cds:
                        pass
                    else:
                        tstart = aa_start*3+offset
                        tend = aa_end*3+offset
                        prot_acc = "%s_%d_%d" % (transcript_id, tstart, tend)
                        tbed = bed.trim(tstart, tend)
                        if args.all or unique_prot(tbed, prot):
                            translate_count += 1
                            tbed.name = prot_acc
                            write_translation(tbed, bed.name, prot)
                    aa_start = aa_end + 1
        return translate_count

    if input_rdr:
        translation_count = 0
        transcript_count = 0
        for i, bedline in enumerate(input_rdr):
            try:
                bed = bed_from_line(bedline, ensembl=args.ensembl,
                                    seq_column=args.column)
                if bed is None:
                    continue
                transcript_count += 1
                if bed.biotype and biotypea and bed.biotype not in biotypea:
                    continue
                if filter_by_regions(bed):
                    translation_count += translate_bed(bed)
            except Exception as e:
                print("BED format Error: line %d: %s\n%s"
                      % (i, bedline, e), file=sys.stderr)
                break
        if args.debug or args.verbose:
            print("transcripts: %d\ttranslations: %d"
                  % (transcript_count, translation_count), file=sys.stderr)


if __name__ == "__main__":
    __main__()
