import os, argparse, shutil, time
from ConfigParser import ConfigParser, NoOptionError

"""Simulate observations using CASA
"""

# Logger class
class Logger:

    def __init__(self, file_name=None):
        self.log = casalog
        if file_name:
            self.log.setlogfile(file_name)

    def info(self, msg, *args, **kwargs):
        self.log.post(msg % args, 'INFO', kwargs.get('name', __name__))

    def warn(self, msg, *args, **kwargs):
        self.log.post(msg % args, 'WARN', kwargs.get('name', __name__))

    def error(self, msg, *args, **kwargs):
        self.log.post(msg % args, 'SEVERE', kwargs.get('name', __name__))

#logger = Logger(__file__.replace('.py', '.log'))
logger = Logger()

# argparse action
class ValidatePath(argparse.Action):

    def __call__(self, parser, namespace, values, option_string=None):
        path = os.path.realpath(os.path.expanduser(os.path.expandvars(values)))
        if os.path.isfile(path):
            setattr(namespace, self.dest, path)
        elif os.path.isdir(path):
            setattr(namespace, self.dest, path)
        else:
            dirname, filename = os.path.split(values)
            if os.path.isdir(dirname):
                pass
            else:
                logger.warn('Creating dir tree:s', dirname, name=__name__)
                os.makedirs(dirname)
            logger.warn('File %s does not exist', filename, name=__name__)
            setattr(namespace, self.dest, path)


def simulate(conf, image, config):
    """Simulate the observation.

    Parameters:
        conf (str): name of the configuration.
        image (str): name of the raw image.
        config (ConfigParser): configuration of the observations.
    """
    try:
        setpointing = config.getboolean(conf,'setpointing')
    except NoOptionError:
        setpointing = True
    try:
        indirection = str(config.get(conf,'indirection'))
    except NoOptionError:
        indirection = ''
    try:
        incenter = str(config.get(conf,'incenter'))
    except NoOptionError:
        incenter = ''
    try:
        inwidth = str(config.get(conf,'inwidth'))
    except NoOptionError:
        inwidth = ''
    try:
        thermalnoise = str(config.get(conf,'thermalnoise'))
    except NoOptionError:
        thermalnoise = ''
    ptgfile = str(config.get(conf,'ptgfile'))
    assert os.path.isfile(ptgfile)
    antennae = os.getenv("CASAPATH").split()[0]+"/data/alma/simmos/"
    antennalist = str(config.get(conf,'antennalist'))
    if antennalist not in antennae:
        assert os.path.isfile(antennalist)

    if setpointing:
        simobserve(project=conf, skymodel=image, indirection=indirection,
                incenter=incenter, inwidth=inwidth, setpointings=setpointing,
                ptgfile=ptgfile, integration=str(config.get(conf,'integration')),
                maptype=str(config.get(conf,'maptype')),
                calflux=str(config.get(conf,'calflux')), 
                refdate=str(config.get(conf,'refdate')),
                totaltime=str(config.get(conf,'totaltime')),
                direction=str(config.get(conf,'indirection')),
                antennalist=antennalist, 
                sdantlist='', sdant=config.getint(conf,'sdant'),
                thermalnoise=thermalnoise,
                user_pwv=config.getfloat(conf,'user_pwv'), 
                t_sky=config.getfloat(conf,'t_sky'),
                t_ground=config.getfloat(conf,'t_ground'),
                tau0=config.getfloat(conf,'tau0'), 
                verbose=config.getboolean(conf,'verbose'), 
                mapsize=['',''], hourangle=str(config.get(conf,'hourangle')),
                graphics=str(config.get(conf,'graphics')))
    else:
        simobserve(project=conf, skymodel=image, indirection=indirection,
                incenter=incenter, inwidth=inwidth, setpointings=setpointing,
                ptgfile=ptgfile, integration=str(config.get(conf,'integration')),
                calflux=str(config.get(conf,'calflux')), 
                refdate=str(config.get(conf,'refdate')), 
                totaltime=str(config.get(conf,'totaltime')),
                antennalist=antennalist, 
                sdantlist='', sdant=config.getint(conf,'sdant'),
                thermalnoise=thermalnoise,
                user_pwv=config.getfloat(conf,'user_pwv'), 
                t_sky=config.getfloat(conf,'t_sky'),
                t_ground=config.getfloat(conf,'t_ground'),
                tau0=config.getfloat(conf,'tau0'), 
                verbose=config.getboolean(conf,'verbose'), 
                hourangle=str(config.get(conf,'hourangle')),
                graphics=str(config.get(conf,'graphics')))
    try:
        shutil.move('simobserve.last', conf)
    except:
        pass

    if thermalnoise: #config.get(conf,'thermalnoise'):
        return str('%s/%s.%s.noisy.ms' % (conf, conf,
            str(config.get(conf,'antenna'))))
    else:
        return str('%s/%s.%s.ms' % (conf, conf, str(config.get(conf,'antenna'))))

def run_simanalyze(conf, vis, config, interactive=False):
    """Clean the simulated visibilities.

    Parameters:
        conf (str): name of the configuration.
        vis (str): name of the visibility ms.
        config (ConfigParser): observation configuration.
        interactive (bool, default=False): interactive cleaning.
    """
    def clean_dir(project, dirname='./'):
        ls = glob.glob(dirname + '%s/%s.concat.*' % (project, project))
        msname = '%s.concat.ms' % project
        for f in ls:
            if os.path.basename(f) != msname:
                shutil.rmtree(f)

    while True:
        try:
            mask = str(config.get(conf,'mask')).split(',')
        except NoOptionError:
            mask = ''
        if len(mask)!=4:
            mask = []
        else:
            mask = map(int, mask)
        try:
            niter = config.getint(conf,'niter')
        except NoOptionError:
            niter = 50000
        try:
            cell = str(config.get(conf,'cell'))
        except NoOptionError:
            cell = ''
        try:
            threshold = str(config.get(conf,'threshold'))
        except NoOptionError:
            threshold = "5e-4Jy"
        try:
            modelimage = str(config.get(conf,'modelimage'))
        except NoOptionError:
            modelimage = ''
        simanalyze(project=conf, imagename='', vis=vis, niter=niter,
                cell=cell, analyze=True, showuv=False, showpsf=True, 
                showmodel=True, mask=mask, threshold=threshold,
                modelimage=modelimage, imsize=0, interactive=interactive,
                imdirection=str(config.get(conf,'indirection')), 
                weighting=str(config.get(conf,'weighting')),
                showconvolved=True, showclean=True, showresidual=True,
                showdifference=True, showfidelity=True,
                verbose=config.getboolean(conf,'verbose'),
                graphics='both') #str(conf['graphics']))
        if not interactive:
            break
        else:
            ans = raw_input('Converged?: ')
            if ans=='y' or ans=='Y':
                logger.info('Breaking interactive cleaning now', name=__name__)
                break
            else:
                clean_dir(conf)
                logger.info('Running again', name=__name__)
                continue

    try:
        shutil.move('simanalyze.last', conf)
    except:
        pass

def concatenate_vis(mss, conf, weights=None):
    """Concatenate the visibilities.

    Parameters:
        mss (list): list of the ms to concatenate.
        conf (str): name of the configuration.
        weights (list, default=None): weight of each ms.
    """
    if len(mss)>1:
        concatvis = str('%s/%s.concat.ms' % (conf, conf))
        if os.path.isdir(concatvis):
            shutil.rmtree(concatvis)
        concat(vis=mss, concatvis=concatvis, visweightscale=weights)
    else:
        concatvis = mss[0]
    return concatvis

def save_image(project, antenna, filename, img_type='image', restfreq=None):
    """Save an image as a fits file.

    Parameters:
        project (str): name of the project directory.
        antenna (str): antenna configuration to save.
        filename (str): output fits file name.
        img_type (str, default=image): type of image to save.
        restfreq (str, default=None): rest frequency of the observation.
    """
    if restfreq:
        imhead(imagename=project+'/%s.%s.%s' % (project, antenna, img_type), 
               mode="put", hdkey="restfreq", hdvalue=restfreq)
    img = ia.open(project+'/%s.%s.%s' % (project, antenna, img_type))
    #check_dir_tree(os.path.dirname(filename))
    ia.tofits(str(filename), overwrite=True)

def main():
    # Args parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', nargs=1,
            help='Casa parameter')
    parser.add_argument('--nologger', action='store_true',
            help='Casa nologger arg')
    parser.add_argument('--nogui', action='store_true',
            help='Casa nogui arg')
    parser.add_argument('--log2term', action='store_true',
            help='Casa log2term arg')
    parser.add_argument('-i', action='store_true',
            help='Interactive cleaning')
    parser.add_argument('config', type=str, action=ValidatePath,
            help='Configuration file name')
    parser.add_argument('image', type=str, action=ValidatePath,
            help='Input image file name')
    parser.add_argument('output', type=str, action=ValidatePath,
            help='Output image file name')
    args = parser.parse_args()

    #logger = Logger(__file__.replace('.py', '.log'))
    logger = Logger(file_name=os.path.join(os.path.dirname(args.output), 
        'casa_simulator.log'))

    # Load configuration files
    config = ConfigParser()
    config.read(args.config)

    # Simulate for each configuration
    concat = []
    weights = []
    for conf in config.sections():
        logger.info('Simulating observations for section: %s', conf,
                name=__name__)
        concat += [simulate(conf, args.image, config)]
        try:
            weights += [config.getfloat(conf,'visweight')]
        except NoOptionError:
            weights += [1.]

    # Concatenate
    logger.info('Concat: %r', concat, name=__name__)
    conf = config.sections()[-1]
    concatenated = concatenate_vis(concat, conf, weights=weights) 

    # Save fits file
    filename = args.output
    restfreq = str(config.get(conf,'restfreq'))
    run_simanalyze(conf, concatenated, config, interactive=args.i)
    save_image(conf, 'concat', filename, restfreq=restfreq)

    # Clean
    #time.sleep(60)
    #for conf in config.sections():
    #    #shutil.rmtree(conf)
    #    sig = os.system('rm -rf CONF*')
    #    print sig
    print 'Finished'

if __name__=='__main__':
    main()
