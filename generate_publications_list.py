import urllib, urllib.request, json

citation_stats = {'max': 0, 'mean': 0, 'total': 0, 'N': 0}

# Use inspirehep API to get publication information
def get_record(recid: str):
    inspire_result = dict()
    inspire_query = 'https://inspirehep.net/api/literature'
    inspire_query += '/' + recid
    with urllib.request.urlopen(inspire_query) as req:
            inspire_result = json.load(req)
    return inspire_result

# generate a list of publication strings
def generate_publication_lists(publication_id_list: list) -> list:
    global citation_stats
    pub_str_list = []
    thesis_id_list= []
    for pub_id in publication_id_list:
        record = get_record(str(pub_id))

        metadata = record['metadata']

        citations = metadata['citation_count']

        author_Str = ", ".join([author['full_name'] for author in metadata['authors']])

        title = metadata['titles'][0]['title']

        refurl = None
        if metadata['document_type'][0] == 'thesis':
            thesis_id_list.append(pub_id)
            continue
        elif 'dois' in metadata:
            DOI = metadata['dois'][0]['value']
            refurl = rf"[DOI:{DOI}](https://doi.org/{DOI})"
        elif 'arxiv_eprints' in metadata:
            arxiv_id = metadata['arxiv_eprints'][0]['value']
            refurl = rf"[arxiv:{arxiv_id}](https://arxiv.org/abs//{arxiv_id})"
        else:
            print(f"error: {pub_id}")
            continue

        output_str = f"-# ***{title}***; {author_Str}; {refurl}; Cites: {citations}  \n"
        pub_str_list.append(output_str)

        # Update statistics dictionary
        citation_stats['max'] = max(citation_stats['max'], citations)
        citation_stats['total'] += citations
        citation_stats['N'] += 1
    return pub_str_list, thesis_id_list

# generate a list of thesis strings
def generate_thesis_list(publication_id_list: list) -> list:
    thesis_str_list = []
    for pub_id in publication_id_list:
        record = get_record(str(pub_id))

        metadata = record['metadata']
        author_data = metadata['authors'][0]

        author_Str = author_data['full_name']

        affilitation = 'Unknown Affiliation'
        if 'affiliations' in author_data:
            affilitation = author_data['affiliations'][0].get('value', 'Unknown Affiliation')
        elif 'raw_affiliations' in author_data:
            affilitation = author_data['raw_affiliations'][0].get('value', 'Unknown Affiliation')


        title = metadata['titles'][0]['title']

        refurl = f"[InspireHep Link](https://inspirehep.net/literature/{pub_id})"

        output_str = f"-# ***{title}***; {author_Str}; *{affilitation}*; {refurl}   \n"
        thesis_str_list.append(output_str)

    return thesis_str_list

"""
This list is defined manually so as to ensure only publications where FUKA was used to enable scientific discovery.
Therefore, it must be updated periodically.
"""
publication_id_list = [
    3182792,
    3181098,
    3179986,
    3178892,
    3178413,
    3167914,
    3167532,
    3165053,
    3163010,
    3151473,
    3136094,
    3136359,
    3127729,
    3130476,
    3112069,
    3110479,
    3142790,
    3098884,
    3097845,
    3093918,
    3086228,
    3062871,
    2974355,
    2960930,
    2959682,
    2958519,
    2954489,
    2926088,
    2945968,
    2939416,
    2928739,
    2912216,
    2911149,
    2900920,
    2893312,
    2874509,
    2871942,
    2857361,
    2856251,
    2854129,
    2854249,
    2852178,
    2838582,
    2831105,
    2828866,
    2827207,
    2827218,
    2812089,
    2811084,
    2810828,
    2796618,
    2794833,
    2791235,
    2789240,
    2781693,
    2780751,
    2769170,
    2765782,
    2759492,
    2816510,
    2738710,
    2748332,
    2719885,
    2714358,
    2711854,
    2695528,
    2693560,
    2678734,
    2734507,
    2668036,
    2653780,
    2642224,
    2641410,
    2635821,
    2593596,
    2080528,
    2009344,
    1915700,
    1868144,
    1835650,
]

pub_str_list, thesis_id_list = generate_publication_lists(publication_id_list)
thesis_str_list = generate_thesis_list(thesis_id_list)
citation_stats['mean'] = citation_stats['total'] / citation_stats['N']

output_file = "FUKA_ENABLED_PUBLICATIONS.md"

with open(output_file, 'w') as f:
    f.write(f'{r"\page"} fukascience Research Enabled by FUKA\n\n')

    f.write('## Known academic theses enabled by the FUKA\n\n')
    for thesis_str in thesis_str_list:
        f.write(thesis_str)

    f.write('\n\n## Known scientific articles enabled by the FUKA\n\n')
    for pub_str in pub_str_list:
        f.write(pub_str)

    f.write('\n\n### Citation Statistics\n\n')
    f.write(f"\nAverage citations of FUKA enabled publications: {citation_stats['mean']:1.0f}  ")
    f.write(f"\nMaximum citations of FUKA enabled publications: {citation_stats['max']}  ")

    f.write('\n\nThis page has been generated using the [Inspire REST API](https://github.com/inspirehep/rest-api-doc)  \n')

    from datetime import datetime
    f.write(f"Last updated: {datetime.today().strftime('%Y-%m-%d')}")