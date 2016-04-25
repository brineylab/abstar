.. _APIexamples:

API Examples
============

AbStar and AbTools both expose a public API containing many of the core functions.
This makes it reasonably straightforward to build custom pipelines that include
several AbStar/AbTools components or integrate these tools with third-party tools. 
A few simple examples are shown below.

Case #1
-------
Sequencing data consists of an Illumina MiSeq run on human samples, with the raw data
stored in BaseSpace (project ID: 123456789). Samples are indexed, so each sample will 
be downloaded from BaseSpace as a separate pair of read files. We'd like to do several things:

  - get a FASTQC report on the raw data
  - remove adapters
  - quality trim
  - get another FASTQC report on the cleaned data
  - merge paired reads
  - annotate with AbStar

::

    import os

    import abstar
    from abstar.utils import basespace, pandaseq

    PROJECT_DIR = '/path/to/project'
    PROJECT_ID = '123456789'

    # download data from BaseSpace
    bs_dir = os.path.join(PROJECT_DIR, 'raw_data')
    basespace.download(bs_dir, project_id=PROJECT_ID)

    # FASTQC on the raw data
    fastqc1_dir = os.path.join(PROJECT_DIR, 'fastqc-pre')
    abstar.fastqc(bs_dir, output=fastqc1_dir)

    # adapter trimming
    adapter_dir = os.path.join(PROJECT_DIR, 'adapter_trimed')
    adapters = '/path/to/adapters.fasta'
    abstar.adapter_trim(bs_dir, output=adapter_dir, adapter_both=adapters)

    # quality trimming
    quality_dir = os.path.join(PROJECT_DIR, 'quality_trimed')
    abstar.quality_trim(adapter_dir, output=quality_dir)

    # FASTQC on the cleaned data
    fastqc2_dir = os.path.join(PROJECT_DIR, 'fastqc-post')
    abstar.fastqc(quality_dir, output=fastqc2_dir)

    # read merging
    merged_dir = os.path.join(PROJECT_DIR, 'merged')
    pandaseq.run(quality_dir, merged_dir)

    # run AbStar
    temp_dir = os.path.join(PROJECT_DIR, 'temp')
    json_dir = os.path.join(PROJECT_DIR, 'json')
    abstar.run(input=merged_dir,
               temp=temp_dir,
               output=json_dir)



Case #2
-------
Sequencing data is a directory of single-read FASTQ files that have already been quality/adapter trimmed. 
We'd like to do the following:

  - get a FASTQC report
  - annotate with AbStar
  - import the JSONs into a MongoDB database named ``MyDatabase``

Our FASTQ file names are formatted as: ``SampleNumber-SampleName.fastq``, which means the AbStar output
file name would be ``SampleNumber-SampleName.json``. We'd like the corresponding MongoDB collection 
to just be named ``SampleName``.

::

    import os

    import abstar
    from abstar.utils import mongoimport

    PROJECT_DIR = '/path/to/project'
    FASTQ_DIR = '/path/to/fastqs'

    MONGO_IP = '123.45.67.89'
    MONGO_PORT = 27017
    MONGO_USER = 'MyUsername'
    MONGO_PASS = 'Secr3t'

    # FASTQC on the input data
    fastqc_dir = os.path.join(PROJECT_DIR, 'fastqc')
    abstar.fastqc(FASTQ_DIR, output=fastqc_dir)

    # run AbStar
    temp_dir = os.path.join(PROJECT_DIR, 'temp')
    json_dir = os.path.join(PROJECT_DIR, 'json')
    abstar.run(input=FASTQ_DIR,
               temp=temp_dir,
               output=json_dir)

    # import into MongoDB
    mongoimport.run(ip=MONGO_IP,
                    port=MONGO_PORT
                    user=MONGO_USER,
                    password=MONGO_PASS,
                    input=json_dir,
                    db='MyDatabase'
                    delim1='-',
                    delim2='.')


Case #3
-------
Now we'd like to use AbStar as part of an analysis script in which sequence annotation 
isn't the primary output. In the previous
examples, we started with raw(ish) sequence data and ended with either a directory of 
JSON files or a MongoDB database populated with AbStar output. In this case, we're 
going to start with a MongoDB database, query that database for some sequences, and 
generate the unmutated common ancestor (UCA). We'd like to annotate the UCA sequence 
inline (as part of the script) so that we can do world-changing things with the 
annotated UCA later in our script. For simplicity's sake, we're querying a local MongoDB 
database that doesn't have authentication enabled, although ``abtools.mongodb`` can 
work with remote MongoDB servers that require authentication.

::

    import abstar

    from abtools import mongodb
    from abtools.sequence import Sequence

    DB_NAME = 'MyDatabase'
    COLLECTION_NAME = 'MyCollection'

    def get_sequences(db_name, collection_name):
        db = mongodb.get_db(db_name)
        c = db[collection]
        seqs = c.find({'chain': 'heavy'})
        return [Sequence(s) for s in seqs]

    def calculate_uca(sequences):
        #
        # code to calculate the UCA sequence, as a string
        #
        return uca

    # get sequences, calculate the UCA
    sequences = get_sequences(DB_NAME, COLLECTION_NAME)
    uca_seq = calculate_uca(sequences)

    # run AbStar on the UCA, returns an AbTools Sequence object
    uca = abstar.run(['UCA', uca_seq])

    # do amazing, world-changing things with the UCA
    # ...
    # ...
    # ... 
