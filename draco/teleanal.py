from astroquery.astrometry_net import AstrometryNet
from astroquery.exceptions  import TimeoutError
from astropy.wcs import WCS

ast = AstrometryNet()
ast.api_key = 'jefjaaeekvrszphp'

__all__ = ['plate_solve']

def  plate_solve(data, image_file_path, write_fits = False, write_name = ''):

        """
        plate_solve takes fits image file data and the corresponding file string to the data and calls nova.astrometry.net to obtain and then integrate WCS into the HDU.

        Parameters
        ----------
        data: CCDdata
                CCDdata of the image to be solved.
        image_file_path: string
                Path to the image to be solved.
        write_fits: Boolean
                Whether to save the resulting data.
        write_name: string
                File name to save the resulting data under.
        Returns
        -------
        wcs_hdu: CCDdata
                The CCDdata file with a wcs header.
        """

        try_again = True
        submission_id = None
        num = 0

        while try_again:
                try:
                        if not submission_id:
                                wcs_header = ast.solve_from_image(image_file_path, force_image_upload=True, submission_id=submission_id, solve_timeout=300)
                        else:
                                print('Monitoring: try #', num)
                                wcs_header = ast.monitor_submission(submission_id, solve_timeout=120)
                except TimeoutError as e:
                        print(TimeoutError)
                        num = num + 1
                        print('Timed out: try #', num)
                        submission_id = e.args[1]


                if wcs_header != None:
                        # got a result, so terminate while loop
                        try_again = False
        if wcs_header:
                # Code to execute when solve succeeds
                print('Solve succeeded! :)')
                print(wcs_header.tostring( sep='\n', endcard=True, padding=True))  # items, keys, values
                wcs_hdu = data
                wcs_hdu.header = wcs_header
                if write_fits:
                        if write_name != '':
                                wcs_hdu.write(write_name + '.fits')
                        else:
                                print('The name to write this file was not provided. File saved under \'draco_wcs.fits\'')
                                wcs_hdu.write('draco_wcs.fits', overwrite = True)
                return wcs_hdu
        else:
                # Code to execute when solve fails
                print('Solve failed! :(')
                return 
