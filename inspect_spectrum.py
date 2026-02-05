#!/usr/bin/env python
"""
Spectrum inspection tool for NIRPS and SPIROU data.

Usage:
    python inspect_spectrum.py                      # All orders PDF
    python inspect_spectrum.py --order 50           # Single order
    python inspect_spectrum.py --order 40 60        # Orders 40-60
"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.table import Table
from scipy.interpolate import InterpolatedUnivariateSpline as IUS
import os
import argparse
from datetime import datetime


def doppler(wave, velocity):
    """Apply Doppler shift to wavelength array."""
    c = 299792.458  # speed of light in km/s
    return wave * np.sqrt((1 - velocity / c) / (1 + velocity / c))


def robust_polyfit(x, y, degree=4, sigma_clip=3.0, max_iter=10):
    """
    Perform robust polynomial fit with iterative sigma clipping.
    
    Parameters
    ----------
    x : array
        Independent variable
    y : array
        Dependent variable
    degree : int
        Polynomial degree (default 4)
    sigma_clip : float
        Number of sigma for outlier rejection (default 3.0)
    max_iter : int
        Maximum number of iterations (default 10)
    
    Returns
    -------
    fit_values : array
        Polynomial evaluated at x positions
    """
    # Start with finite values only
    mask = np.isfinite(x) & np.isfinite(y)
    
    for iteration in range(max_iter):
        if np.sum(mask) < degree + 1:
            # Not enough points for fit
            return np.ones_like(y) * np.nan
        
        # Fit polynomial to current valid points
        coeffs = np.polyfit(x[mask], y[mask], degree)
        fit = np.polyval(coeffs, x)
        
        # Compute residuals and sigma
        residuals = y - fit
        sigma = np.nanstd(residuals[mask])
        
        if sigma == 0:
            break
        
        # Update mask - reject outliers
        new_mask = np.isfinite(x) & np.isfinite(y) & (np.abs(residuals) < sigma_clip * sigma)
        
        # Check for convergence
        if np.array_equal(mask, new_mask):
            break
        
        mask = new_mask
    
    return fit


def get_file_type(fits_file):
    """
    Determine if the file is an APERO or ESO file.
    
    Parameters
    ----------
    fits_file : str
        Path to the FITS file
    
    Returns
    -------
    str
        'APERO' if file ends with 't.fits'
        'ESO' if file basename starts with 'r.'
        None if file type cannot be determined
    """
    basename = os.path.basename(fits_file)
    
    if fits_file.endswith('t.fits'):
        return 'APERO'
    elif basename.startswith('r.'):
        return 'ESO'
    else:
        return None


def get_oh_lines(wave_min, wave_max):
    """Download and return OH emission lines in the specified range."""
    # Store OH data in reference_data folder
    os.makedirs('reference_data', exist_ok=True)
    oh_file = 'reference_data/tablea1.dat'
    
    url = 'https://cdsarc.cds.unistra.fr/ftp/J/A+A/581/A47/table1.dat'
    if not os.path.exists(oh_file):
        import urllib.request
        urllib.request.urlretrieve(url, oh_file)
    
    oh_lines = Table.read(oh_file, format='ascii')
    wave_oh = np.concatenate([
        np.array(oh_lines['col1'].data),
        np.array(oh_lines['col3'].data)
    ]) / 10.0  # convert Angstrom to nm
    label_oh = np.concatenate([oh_lines['col2'].data, oh_lines['col4'].data])
    
    keep = (wave_oh >= wave_min) & (wave_oh <= wave_max)
    return wave_oh[keep], label_oh[keep]


def load_tapas(tapas_file):
    """Load TAPAS atmospheric transmission data."""
    if not os.path.exists(tapas_file):
        # Try to download from default URL
        # Ensure reference_data folder exists
        os.makedirs('reference_data', exist_ok=True)
        
        url = 'http://206.12.93.77/ari/data/lbl/tapas/tapas_lbl.fits'
        print(f"TAPAS file '{tapas_file}' not found locally.")
        print(f"Attempting to download from: {url}")
        
        try:
            import urllib.request
            urllib.request.urlretrieve(url, tapas_file)
            
            # Verify the file was actually downloaded and is valid
            if os.path.exists(tapas_file):
                file_size = os.path.getsize(tapas_file)
                if file_size < 1000:  # Less than 1KB is suspicious
                    os.remove(tapas_file)
                    raise ValueError(f"Downloaded file is too small ({file_size} bytes). Download may have failed.")
                print(f"✓ Successfully downloaded tapas_lbl.fits ({file_size / (1024*1024):.1f} MB)")
            else:
                raise FileNotFoundError("Download completed but file not found.")
                
        except Exception as e:
            raise FileNotFoundError(
                f"Failed to download TAPAS file from {url}\n"
                f"Error: {e}\n\n"
                f"The online storage may be temporarily unavailable.\n"
                f"Please try:\n"
                f"  1. Copy 'tapas_lbl.fits' from your LBL 'models/' folder to this directory\n"
                f"  2. Contact the maintainer if the problem persists"
            )
    
    return Table.read(tapas_file)


def load_template(template_file):
    """Load template spectrum and return interpolator."""
    tbl_template = Table(fits.getdata(template_file))
    wave_template = tbl_template['wavelength'].data
    flux_template = tbl_template['flux'].data
    
    keep = np.isfinite(wave_template) & np.isfinite(flux_template)

    return IUS(wave_template[keep], flux_template[keep], k=1),IUS(wave_template,np.array(keep,dtype = float),k=1)


def load_spectrum(fits_file, instrument='NIRPS'):
    """Load spectrum data from FITS file."""
    if instrument == 'NIRPS':
        fiber_setup = 'A'
    elif instrument == 'SPIROU':
        fiber_setup = 'AB'
    else:
        raise ValueError("Instrument must be either 'NIRPS' or 'SPIROU'.")
    
    file_type = get_file_type(fits_file)
    
    # ESO files have different header structure - BERV is in primary header
    if file_type == 'ESO':
        hdr_primary = fits.getheader(fits_file, ext=0)
        hdr = fits.getheader(fits_file, ext=1)
        berv = hdr_primary['ESO QC BERV']
        # Get systemic velocity from primary header
        try:
            syst_vel = hdr_primary['ESO TEL TARG RADVEL']
        except KeyError:
            syst_vel = np.nan
    else:
        hdr = fits.getheader(fits_file, ext=1)
        berv = hdr['BERV']
        # Get systemic velocity (target radial velocity)
        try:
            syst_vel = hdr['ESO TEL TARG RADVEL']
        except KeyError:
            syst_vel = np.nan
    
    # Load data - different extensions for APERO vs ESO
    if file_type == 'ESO':
        # ESO files: SCIDATA for flux, WAVEDATA_VAC_BARY for wavelength
        # Note: ESO wavelength is already BERV-corrected (rest frame)
        flux = fits.getdata(fits_file, 'SCIDATA')
        wave = fits.getdata(fits_file, 'WAVEDATA_VAC_BARY')
        # ESO wavelength is in Angstroms, convert to nm
        wave = wave / 10.0
        # ESO files don't have blaze extension - create unity blaze
        blaze = np.ones_like(flux)
        # Flag that wavelength is already in rest frame
        wave_in_rest_frame = True
    else:
        # APERO files: FluxA/FluxAB, WaveA/WaveAB, BlazeA/BlazeAB
        flux = fits.getdata(fits_file, f'Flux{fiber_setup}')
        blaze = fits.getdata(fits_file, f'Blaze{fiber_setup}')
        wave = fits.getdata(fits_file, f'Wave{fiber_setup}')
        wave_in_rest_frame = False
    
    # Convert zeros to NaN - there are no real zeros in spectral data
    flux[flux == 0] = np.nan
    
    # Normalize blaze
    for iord in range(flux.shape[0]):
        blaze[iord] /= np.nanpercentile(blaze[iord], 95)
    
    # Extract target name and date from header
    # For ESO, use primary header; for APERO, use ext=1 header
    if file_type == 'ESO':
        target = hdr_primary.get('OBJECT', 'Unknown')
        date_obs = hdr_primary.get('DATE-OBS', datetime.now().strftime('%Y-%m-%d'))
    else:
        target = hdr.get('OBJECT', 'Unknown')
        date_obs = hdr.get('DATE-OBS', datetime.now().strftime('%Y-%m-%d'))
    
    return {
        'flux': flux,
        'blaze': blaze,
        'wave': wave,
        'berv': berv,
        'syst_vel': syst_vel,
        'v_tot': berv - syst_vel if np.isfinite(syst_vel) else np.nan,
        'target': target,
        'date_obs': date_obs,
        'instrument': instrument,
        'n_orders': flux.shape[0],
        'fits_file': os.path.basename(fits_file),
        'file_type': get_file_type(fits_file),
        'wave_in_rest_frame': wave_in_rest_frame,
        'needs_blaze_correction': (file_type == 'ESO' and '_BLAZE_' in os.path.basename(fits_file)),
    }


def plot_order(spectrum_data, template_interp, order, tapas=None, 
               wave_oh=None, label_oh=None, show_oh=True, show_tapas=True):
    """
    Plot a single spectral order.
    
    Parameters
    ----------
    spectrum_data : dict
        Dictionary with flux, blaze, wave, berv keys
    template_interp : callable
        Interpolator for template spectrum
    order : int
        Order number to plot
    tapas : Table, optional
        TAPAS atmospheric data (full table, will be trimmed)
    wave_oh, label_oh : array, optional
        OH line positions and labels (full arrays, will be trimmed)
    show_oh : bool
        Whether to show OH lines
    show_tapas : bool
        Whether to show TAPAS transmission
    
    Returns
    -------
    fig, ax : matplotlib figure and axes, or (None, None) if order is invalid
    """
    wave = spectrum_data['wave']
    flux = spectrum_data['flux']
    blaze = spectrum_data['blaze']
    berv = spectrum_data['berv']
    wave_in_rest_frame = spectrum_data.get('wave_in_rest_frame', False)
    needs_blaze_correction = spectrum_data.get('needs_blaze_correction', False)
    
    # Get wavelength range for this order
    wave_ord = wave[order]
    
    # Check if flux is all NaN - skip this order
    if np.all(~np.isfinite(flux[order])):
        return None, None
    # For APERO files: divide flux by blaze to remove blaze shape
    # For ESO files: blaze is unity, no effect
    flux_blaze_corrected = flux[order] / blaze[order]
    
    flux_median = np.nanmedian(flux_blaze_corrected)
    if not np.isfinite(flux_median) or flux_median == 0:
        return None, None
    
    # Create figure
    fig, ax = plt.subplots(figsize=(18, 8), nrows=2, sharex=True)
    
    # Normalize the blaze-corrected flux
    plot_tmp = flux_blaze_corrected / flux_median
    
    # For ESO files, wavelength is already in rest frame (BERV-corrected)
    # For APERO files, we need to Doppler shift to rest frame
    if wave_in_rest_frame:
        template_tmp = template_interp[0](wave_ord)
        template_mask = template_interp[1](wave_ord) > 0.5
    else:
        template_tmp = template_interp[0](doppler(wave_ord, -berv))
        template_mask = template_interp[1](doppler(wave_ord, -berv)) > 0.5
    template_tmp[~template_mask] = np.nan

    # Template is already flat - just normalize by median
    template_tmp /= np.nanmedian(template_tmp)
    
    # Store original flux for background plot (before blaze correction)
    plot_tmp_original = None
    flux_label = 'Normalized Flux'
    
    # Apply blaze correction for ESO files with _BLAZE_ in name
    if needs_blaze_correction:
        # Compute ratio of observed spectrum to template
        ratio = plot_tmp / template_tmp
        
        # Robust polynomial fit to the ratio
        blaze_fit = robust_polyfit(wave_ord, ratio, degree=4, sigma_clip=3.0, max_iter=10)
        
        # Check if fit succeeded
        if np.any(np.isfinite(blaze_fit)):
            # Store original for background plot
            plot_tmp_original = plot_tmp.copy()
            
            # Divide observed spectrum by the fit to match template shape
            plot_tmp = plot_tmp / blaze_fit
            
            # Re-normalize after correction
            plot_tmp /= np.nanmedian(plot_tmp[np.isfinite(plot_tmp)])
            
            # Update label to indicate correction was applied
            flux_label = 'Normalized Flux (blaze-corrected)'
    
    # Find valid domain based on where the normalized flux is finite
    # Use indices to find first and last valid points (trimming NaN edges)
    flux_valid = np.isfinite(plot_tmp)
    if not np.any(flux_valid):
        plt.close(fig)
        return None, None
    
    valid_indices = np.where(flux_valid)[0]
    idx_min, idx_max = valid_indices[0], valid_indices[-1]
    wave_min, wave_max = wave_ord[idx_min], wave_ord[idx_max]
    
    # Plot original (non-corrected) spectrum in background if blaze correction was applied
    if plot_tmp_original is not None:
        ax[0].plot(wave_ord, plot_tmp_original, alpha=0.3, color='grey', 
                   label='Original (uncorrected)', zorder=1, rasterized=True)
    
    ax[0].plot(wave_ord, plot_tmp, alpha=0.7, color='black', label=flux_label, zorder=2, rasterized=True)
    ax[0].plot(wave_ord, template_tmp, alpha=0.7, color='red', label='Template', zorder=3, rasterized=True)
    ax[1].plot(wave_ord, plot_tmp - template_tmp, color='orange', label='Residuals', rasterized=True)
    
    # Add OH lines (trim to order range)
    if show_oh and wave_oh is not None:
        keep = (wave_oh >= wave_min) & (wave_oh <= wave_max)
        for w, lab in zip(wave_oh[keep], label_oh[keep]):
            ax[0].axvline(w, color='cyan', alpha=1.0, linewidth=1.0)
            ax[1].axvline(w, color='cyan', alpha=1.0, linewidth=1.0)
            ax[0].text(w, 0.05, lab, rotation=90, verticalalignment='bottom',
                      horizontalalignment='right', fontsize=5, alpha=0.5)
    
    # Add TAPAS (trim to order range)
    if show_tapas and tapas is not None:
        keep = (tapas['WAVELENGTH'] >= wave_min) & (tapas['WAVELENGTH'] <= wave_max)
        tapas_trim = tapas[keep]
        if len(tapas_trim) > 0:
            ax[0].plot(tapas_trim['WAVELENGTH'], tapas_trim['ABSO_OTHERS'], 
                      color='green', alpha=0.9, label='TAPAS Other Gases', linewidth=1.0, rasterized=True)
            ax[0].plot(tapas_trim['WAVELENGTH'], tapas_trim['ABSO_WATER'],
                      color='blue', alpha=0.9, label='TAPAS Water', linewidth=1.0, rasterized=True)
    
    ax[0].legend(loc='upper right', fontsize=8)
    ax[0].set_xlim(wave_min, wave_max)
    
    # Set y-axis limits for flux plot
    if plot_tmp_original is not None:
        # For blaze-corrected ESO data, use 1.2x the 90th percentile of original
        ymax_flux = 1.2 * np.nanpercentile(plot_tmp_original, 90)
        ax[0].set_ylim(0, ymax_flux)
    else:
        ax[0].set_ylim(bottom=0)
    
    ax[1].set_xlabel('Wavelength (nm)')
    ax[0].set_ylabel('Normalized Flux')
    ax[1].set_ylabel('Residuals')
    
    # Set residual ylim - symmetric and capped at [-0.5, 0.5]
    residuals = plot_tmp - template_tmp
    res_valid = residuals[np.isfinite(residuals)]
    if len(res_valid) > 0:
        ymin, ymax = np.min(res_valid), np.max(res_valid)
        # Cap at ±0.5
        ymin = max(ymin, -0.5)
        ymax = min(ymax, 0.5)
        # Force symmetry using mean of absolute values
        ylim = (np.abs(ymin) + np.abs(ymax)) / 2
        ax[1].set_ylim(-ylim, ylim)
    else:
        ax[1].set_ylim(-0.5, 0.5)
    
    # Title
    title = f"{spectrum_data['instrument']} | {spectrum_data['target']} | {spectrum_data['date_obs']} | Order {order} | {wave_min:.1f}-{wave_max:.1f} nm"
    ax[0].set_title(title)
    
    ax[0].grid(alpha=0.3)
    ax[1].grid(alpha=0.3)
    
    return fig, ax


def generate_pdf_name(spectrum_data, order_min=None, order_max=None):
    """Generate smart PDF filename."""
    # Format timestamp - replace colons with dashes for filename safety
    timestamp = spectrum_data['date_obs'].replace(':', '-')
    
    parts = [
        spectrum_data['instrument'],
        spectrum_data['file_type'],  # ESO or APERO
        spectrum_data['target'].replace(' ', '_'),
        timestamp,
    ]
    
    if order_min is not None and order_max is not None:
        if order_min == order_max:
            parts.append(f'order{order_min}')
        else:
            parts.append(f'orders{order_min}-{order_max}')
    else:
        parts.append('all_orders')
    
    return '_'.join(parts) + '.pdf'


def create_summary_page(spectrum_data, template_file, fits_file, order_min, order_max):
    """
    Create a summary page with observation information.
    
    Parameters
    ----------
    spectrum_data : dict
        Dictionary with spectrum metadata
    template_file : str
        Path to template file
    fits_file : str
        Path to input FITS file
    order_min, order_max : int
        Range of orders being plotted
    
    Returns
    -------
    fig : matplotlib figure
    """
    fig, ax = plt.subplots(figsize=(18, 8))
    ax.axis('off')
    
    # Title
    ax.text(0.5, 0.95, f"Spectrum Inspection Report", 
            fontsize=24, fontweight='bold', ha='center', va='top', transform=ax.transAxes)
    ax.text(0.5, 0.88, f"{spectrum_data['instrument']} - {spectrum_data['target']}",
            fontsize=18, ha='center', va='top', transform=ax.transAxes)
    
    # Left column - Observation & Files
    left_lines = [
        "OBSERVATION DETAILS",
        "=" * 35,
        "",
        f"Target:        {spectrum_data['target']}",
        f"Instrument:    {spectrum_data['instrument']}",
        f"Date:          {spectrum_data['date_obs']}",
        "",
        "",
        "FILES",
        "=" * 35,
        "",
        f"Spectrum:",
        f"  {os.path.basename(fits_file)}",
        "",
        f"Template:",
        f"  {os.path.basename(template_file)}",
    ]
    
    # Right column - Velocities & Orders
    syst_vel_str = f"{spectrum_data['syst_vel']:.3f} km/s" if np.isfinite(spectrum_data['syst_vel']) else "N/A"
    v_tot_str = f"{spectrum_data['v_tot']:.3f} km/s" if np.isfinite(spectrum_data['v_tot']) else "N/A"
    
    right_lines = [
        "VELOCITIES",
        "=" * 35,
        "",
        f"BERV:          {spectrum_data['berv']:.3f} km/s",
        f"Systemic (Vsys): {syst_vel_str}",
        f"V_tot (BERV-Vsys): {v_tot_str}",
        "",
        "",
        "SPECTRAL ORDERS",
        "=" * 35,
        "",
        f"Total orders:  {spectrum_data['n_orders']}",
        f"In this PDF:   {order_min} to {order_max}",
        f"               ({order_max - order_min + 1} orders)",
        "",
        f"Wavelength:    {np.nanmin(spectrum_data['wave']):.1f} - {np.nanmax(spectrum_data['wave']):.1f} nm",
    ]
    
    left_text = '\n'.join(left_lines)
    right_text = '\n'.join(right_lines)
    
    # Place columns
    ax.text(0.25, 0.78, left_text, fontsize=12, fontfamily='monospace',
            ha='center', va='top', transform=ax.transAxes)
    ax.text(0.75, 0.78, right_text, fontsize=12, fontfamily='monospace',
            ha='center', va='top', transform=ax.transAxes)
    
    # Footer
    ax.text(0.5, 0.05, f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            fontsize=10, ha='center', va='bottom', transform=ax.transAxes, alpha=0.5)
    
    return fig


def inspect_all_orders(fits_file, template_file, tapas_file='reference_data/tapas_lbl.fits',
                       instrument='NIRPS', order_min=None, order_max=None,
                       show_oh=True, show_tapas=True):
    """
    Create multi-page PDF with one page per spectral order.
    
    Parameters
    ----------
    fits_file : str
        Path to NIRPS/SPIROU FITS file
    template_file : str
        Path to template spectrum
    tapas_file : str
        Path to TAPAS file
    instrument : str
        'NIRPS' or 'SPIROU'
    order_min, order_max : int, optional
        Range of orders to plot. If None, plots all orders.
    show_oh : bool
        Show OH emission lines
    show_tapas : bool
        Show TAPAS transmission
    
    Returns
    -------
    pdf_path : str
        Path to saved PDF
    """
    # Load spectrum data once
    spectrum_data = load_spectrum(fits_file, instrument)
    template_interp = load_template(template_file)
    n_orders = spectrum_data['n_orders']
    
    # Determine order range
    if order_min is None:
        order_min = 0
    if order_max is None:
        order_max = n_orders - 1
    
    # Clamp to valid range
    order_min = max(0, order_min)
    order_max = min(n_orders - 1, order_max)
    
    # Load full OH lines and TAPAS (will trim per order)
    wave_oh, label_oh = None, None
    tapas = None
    
    if show_oh:
        # Get full wavelength range
        full_wave_min = np.nanmin(spectrum_data['wave'])
        full_wave_max = np.nanmax(spectrum_data['wave'])
        wave_oh, label_oh = get_oh_lines(full_wave_min, full_wave_max)
    
    if show_tapas:
        tapas = load_tapas(tapas_file)
    
    # Generate PDF name
    pdf_path = generate_pdf_name(spectrum_data, order_min, order_max)
    
    print(f"Creating PDF with orders {order_min} to {order_max} ({order_max - order_min + 1} orders + summary page)...")
    
    with PdfPages(pdf_path) as pdf:
        # First page: Summary
        fig_summary = create_summary_page(spectrum_data, template_file, fits_file, order_min, order_max)
        pdf.savefig(fig_summary, dpi=150)
        plt.close(fig_summary)
        
        # Order pages
        for order in range(order_min, order_max + 1):
            print(f"  Processing order {order}/{order_max}...", end='\r')
            
            fig, ax = plot_order(
                spectrum_data, template_interp, order,
                tapas=tapas, wave_oh=wave_oh, label_oh=label_oh,
                show_oh=show_oh, show_tapas=show_tapas
            )
            
            # Skip invalid orders
            if fig is None:
                continue
            
            plt.tight_layout()
            pdf.savefig(fig, dpi=150)
            plt.close(fig)
    
    print(f"\nSaved: {pdf_path}")
    return pdf_path


def inspect_single_order(fits_file, template_file, order, tapas_file='reference_data/tapas_lbl.fits',
                         instrument='NIRPS', show_oh=True, show_tapas=True,
                         save_pdf=True, show=True):
    """
    Plot and optionally save a single spectral order.
    
    Parameters
    ----------
    fits_file : str
        Path to NIRPS/SPIROU FITS file
    template_file : str
        Path to template spectrum
    order : int
        Order number to plot
    tapas_file : str
        Path to TAPAS file
    instrument : str
        'NIRPS' or 'SPIROU'
    show_oh : bool
        Show OH emission lines
    show_tapas : bool
        Show TAPAS transmission
    save_pdf : bool
        Save plot to PDF
    show : bool
        Display plot
    
    Returns
    -------
    pdf_path : str or None
        Path to saved PDF if save_pdf=True
    """
    # Load data
    spectrum_data = load_spectrum(fits_file, instrument)
    template_interp = load_template(template_file)
    
    # Load auxiliary data
    wave_ord = spectrum_data['wave'][order]
    wave_min, wave_max = np.nanmin(wave_ord), np.nanmax(wave_ord)
    
    wave_oh, label_oh = get_oh_lines(wave_min, wave_max) if show_oh else (None, None)
    tapas = load_tapas(tapas_file) if show_tapas else None
    
    # Plot
    fig, ax = plot_order(
        spectrum_data, template_interp, order,
        tapas=tapas, wave_oh=wave_oh, label_oh=label_oh,
        show_oh=show_oh, show_tapas=show_tapas
    )
    
    # Check if order is valid
    if fig is None:
        print(f"Order {order} has no valid data, skipping.")
        return None
    
    plt.tight_layout()
    
    # Save PDF
    pdf_path = None
    if save_pdf:
        pdf_path = generate_pdf_name(spectrum_data, order, order)
        fig.savefig(pdf_path, dpi=150)
        print(f"Saved: {pdf_path}")
    
    if show:
        plt.show()
    else:
        plt.close(fig)
    
    return pdf_path


def find_template(fits_file, instrument='NIRPS'):
    """
    Auto-detect template file based on object name from FITS header.
    
    Parameters
    ----------
    fits_file : str
        Path to the spectrum FITS file
    instrument : str
        'NIRPS' or 'SPIROU'
    
    Returns
    -------
    template_path : str
        Path to the template file
    
    Raises
    ------
    FileNotFoundError
        If no matching template is found
    """
    hdr = fits.getheader(fits_file, ext=1)
    obj_name = hdr.get('DRSOBJN', 'Unknown').replace(' ', '_')
    
    if instrument == 'NIRPS':
        fiber = 'A'
    else:
        fiber = 'AB'
    
    template_name = f'Template_s1dv_{obj_name}_sc1d_v_file_{fiber}.fits'
    
    if os.path.exists(template_name):
        print(f"Using auto-detected template: {template_name}")
        return template_name
    else:
        raise FileNotFoundError(
            f"Template file not found. Looked for: {template_name}\n"
            f"Please provide the template file as the second argument, or ensure the file exists."
        )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Inspect NIRPS/SPIROU spectra by order',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  python inspect_spectrum.py spectrum.fits template.fits
  python inspect_spectrum.py spectrum.fits                   # Auto-detect template
  python inspect_spectrum.py spectrum.fits template.fits --order 50
  python inspect_spectrum.py spectrum.fits template.fits --order 40 60
"""
    )
    parser.add_argument('file', help='Input FITS file (e.g., NIRPS.2024-09-28T23:54:06.014t.fits)')
    parser.add_argument('template', nargs='?', default=None,
                        help='Template spectrum file (optional, will auto-detect if not provided)')
    parser.add_argument('--tapas', default='reference_data/tapas_lbl.fits',
                        help='TAPAS transmission file')
    parser.add_argument('--instrument', default='NIRPS', choices=['NIRPS', 'SPIROU'],
                        help='Instrument name')
    parser.add_argument('--order', type=int, nargs='+',
                        help='Order(s) to plot. Single value for one order, two values for range.')
    parser.add_argument('--no-oh', action='store_true', help='Hide OH lines')
    parser.add_argument('--no-tapas', action='store_true', help='Hide TAPAS')
    parser.add_argument('--no-save', action='store_true', help='Do not save PDF')
    parser.add_argument('--no-show', action='store_true', help='Do not display plot')
    
    args = parser.parse_args()
    
    # Check input file exists
    if not os.path.exists(args.file):
        print(f"Error: Input file '{args.file}' not found.")
        exit(1)
    
    # Check that the file is a valid APERO or ESO file
    file_type = get_file_type(args.file)
    if file_type is None:
        print(f"""\n{'='*70}
ERROR: Invalid file type!
{'='*70}

You provided: {args.file}

This script works with two types of files:

1. APERO telluric-corrected files (t.fits):
   - File name MUST end with 't.fits'
   - Example: NIRPS.2024-09-28T23:54:06.014t.fits

2. ESO reduced files:
   - File name MUST start with 'r.'
   - Example: r.NIRPS.2023-08-25T01_29_54.429_S2D_BLAZE_TELL_CORR_A.fits

Files that do NOT work:
  - NIRPS.2024-09-28T23:54:06.014e.fits  ← WRONG (extracted, not telluric-corrected)
  - NIRPS.2024-09-28T23:54:06.014o.fits  ← WRONG (not telluric-corrected)

Please use a valid APERO t.fits or ESO r. file.
{'='*70}\n""")
        exit(1)
    
    print(f"Detected file type: {file_type}")
    
    # Auto-detect template if not provided (only for APERO files)
    if args.template is None:
        if file_type == 'ESO':
            print(f"""\n{'='*70}
ERROR: Template file required for ESO files!
{'='*70}

ESO files do not have reliable OBJECT names for template auto-detection.
You must provide a template file explicitly:

  python inspect_spectrum.py {args.file} <template.fits>

{'='*70}\n""")
            exit(1)
        try:
            args.template = find_template(args.file, args.instrument)
        except FileNotFoundError as e:
            print(f"Error: {e}")
            exit(1)
    elif not os.path.exists(args.template):
        print(f"Error: Template file '{args.template}' not found.")
        exit(1)
    
    if args.order is None:
        # All orders
        inspect_all_orders(
            args.file, args.template, args.tapas,
            instrument=args.instrument,
            show_oh=not args.no_oh,
            show_tapas=not args.no_tapas
        )
    elif len(args.order) == 1:
        # Single order
        inspect_single_order(
            args.file, args.template, args.order[0],
            tapas_file=args.tapas,
            instrument=args.instrument,
            show_oh=not args.no_oh,
            show_tapas=not args.no_tapas,
            save_pdf=not args.no_save,
            show=not args.no_show
        )
    else:
        # Range of orders
        inspect_all_orders(
            args.file, args.template, args.tapas,
            instrument=args.instrument,
            order_min=args.order[0],
            order_max=args.order[1],
            show_oh=not args.no_oh,
            show_tapas=not args.no_tapas
        )
