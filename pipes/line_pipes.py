from myutils.decorrators import register_function

def mollie_base_pipe(config, model_file, PA=0.):
    # Load & save FITS file

    # Rotate & save

    return filename


@register_function
def mollie_casa_line_pv(config, model_file, PA=0., vlsr=0.):
    # Run basic pipeline
    filename = mollie_base_pipe(config, model_file, PA=PA)

    # Run casa

    # Get PV maps

    return data


