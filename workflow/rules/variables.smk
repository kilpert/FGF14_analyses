## outdir
try:
    results = config["outdir"]
except:
    config["outdir"] = "results"
    results = config["outdir"]


## verbose
try:
    verbose = config["verbose"]
except:
    verbose = False


## fastq ext
try:
    config["fastq_ext"]
except:
    config["fastq_ext"] = ".fastq.gz"


## sample whitelist
try:
    config["whitelist"]
except:
    config["whitelist"] = None

## sample blacklist
try:
    config["blacklist"]
except:
    config["blacklist"] = None


## orientations
try:
    orientations = config["orientations"]
except:
    config["orientations"] = ['FF', 'RR']
    orientations = config["orientations"]



## sample
try:
    config["sample"] = {}
    for fastq in glob.glob(config["fastq_dir"]+"/**/*"+config["fastq_ext"], recursive=True):
        sample = os.path.basename(fastq).split(".")[0]
        if config["whitelist"]:
            if sample in config["whitelist"]:
                config["sample"][sample] = fastq
        else:
            config["sample"][sample] = fastq

        if config["blacklist"]:
            if sample in config["blacklist"]:
                config["sample"].pop(sample, None) # remove blacklist key

    samples = sorted(list(config["sample"].keys()))
except:
    pass


## ref
try:
    refs = config["data_repo"]["ref"].keys()
except:
    refs = None

