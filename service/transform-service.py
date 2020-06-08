from flask import Flask, request, Response
import json
import os
import math
import logging
from utils import parse_json_stream, entities_to_json

app = Flask(__name__)

logger = None

easting_property = os.environ.get("EASTING_PROPERTY", "easting")
northing_property = os.environ.get("NORTHING_PROPERTY", "northing")
zone_property = os.environ.get("ZONE_PROPERTY", "zone")
zone_default = os.environ.get("ZONE_DEFAULT", "32")
hemi_property = os.environ.get("HEMI_PROPERTY", "hemi")
hemi_default = os.environ.get("HEMI_DEFAULT", "0")
hemi_northern_value = os.environ.get("HEMI_NORTHERN_VALUE", "0")
lat_property = os.environ.get("LATITUDE_PROPERTY", "lat")
long_property = os.environ.get("LONGITUDE_PROPERTY", "long")
include_latlong = os.environ.get("INCLUDE_LAT_LONG", "False")
latlong_property = os.environ.get("LAT_LONG_PROPERTY", "lat_long")

#################################################################
# All the following code to calculate coordinates is copied     #
# and rewritten as python from                                  #
# http://home.hiwaay.net/~taylorc/toolbox/geography/geoutm.html #
#################################################################

pi = 3.14159265358979
# Ellipsoid model constants (actual values here are for WGS84)
sm_a = 6378137.0
sm_b = 6356752.314
sm_EccSquared = 6.69437999013e-03

UTMScaleFactor = 0.9996


def utm_central_meridian(zone):
    """
    Determines the central meridian for the given UTM zone.

    Inputs:
        zone - An integer value designating the UTM zone, range [1,60].

    Returns:
        The central meridian for the given UTM zone, in radians, or zero
        if the UTM zone parameter is outside the range [1,60].
        Range of the central meridian is the radian equivalent of [-177,+177].
    """
    return deg_to_rad(-183.0 + (zone * 6.0))


def deg_to_rad(deg):
    """
    Converts degrees to radians.
    """
    return deg / 180.0 * pi


def rad_to_deg(rad):
    """
    Converts radians to degrees.
    """
    return rad / pi * 180.0


def footpoint_laititude(y):
    """
    Computes the footpoint latitude for use in converting transverse
    Mercator coordinates to ellipsoidal coordinates.

    Reference: Hoffmann-Wellenhof, B., Lichtenegger, H., and Collins, J.,
        GPS: Theory and Practice, 3rd ed.  New York: Springer-Verlag Wien, 1994.

    Inputs:
        y - The UTM northing coordinate, in meters.

    Returns:
        The footpoint latitude, in radians.
    """
    # Precalculate n (Eq. 10.18)
    n = (sm_a - sm_b) / (sm_a + sm_b)

    # Precalculate alpha_ (Eq. 10.22)
    # (Same as alpha in Eq. 10.17)
    alpha_ = ((sm_a + sm_b) / 2.0) * (1 + (math.pow(n, 2.0) / 4) + (math.pow(n, 4.0) / 64))

    # Precalculate y_ (Eq. 10.23)
    y_ = y / alpha_

    # Precalculate beta_ (Eq. 10.22)
    beta_ = (3.0 * n / 2.0) + (-27.0 * math.pow(n, 3.0) / 32.0) + (269.0 * math.pow(n, 5.0) / 512.0)

    # Precalculate gamma_ (Eq. 10.22)
    gamma_ = (21.0 * math.pow(n, 2.0) / 16.0) + (-55.0 * math.pow(n, 4.0) / 32.0)

    # Precalculate delta_ (Eq. 10.22)
    delta_ = (151.0 * math.pow(n, 3.0) / 96.0) + (-417.0 * math.pow(n, 5.0) / 128.0)

    # Precalculate epsilon_ (Eq. 10.22)
    epsilon_ = (1097.0 * math.pow(n, 4.0) / 512.0)

    # Now calculate the sum of the series (Eq. 10.21)
    result = y_ + (beta_ * math.sin(2.0 * y_)) + (gamma_ * math.sin(4.0 * y_)) + \
             (delta_ * math.sin(6.0 * y_)) + (epsilon_ * math.sin(8.0 * y_))

    return result


def map_xy_to_lat_lon(x, y, lambda0):
    """
    Converts x and y coordinates in the Transverse Mercator projection to
    a latitude/longitude pair.  Note that Transverse Mercator is not
    the same as UTM; a scale factor is required to convert between them.

    Reference: Hoffmann-Wellenhof, B., Lichtenegger, H., and Collins, J.,
        GPS: Theory and Practice, 3rd ed.  New York: Springer-Verlag Wien, 1994.

    Inputs:
        x - The easting of the point, in meters.
        y - The northing of the point, in meters.
        lambda0 - Longitude of the central meridian to be used, in radians.

    Outputs:
        philambda - A 2-element containing the latitude and longitude
                   in radians.

    Returns:
        The function does not return a value.

    Remarks:
       The local variables Nf, nuf2, tf, and tf2 serve the same purpose as
       N, nu2, t, and t2 in MapLatLonToXY, but they are computed with respect
       to the footpoint latitude phif.

       x1frac, x2frac, x2poly, x3poly, etc. are to enhance readability and
       to optimize computations.
    """
    # Get the value of phif, the footpoint latitude.
    phif = footpoint_laititude(y)

    # Precalculate ep2
    ep2 = (math.pow(sm_a, 2.0) - math.pow(sm_b, 2.0)) / math.pow(sm_b, 2.0)

    # Precalculate cos (phif)
    cf = math.cos(phif)

    # Precalculate nuf2
    nuf2 = ep2 * math.pow(cf, 2.0)

    # Precalculate Nf and initialize Nfpow
    Nf = math.pow(sm_a, 2.0) / (sm_b * math.sqrt(1 + nuf2))
    Nfpow = Nf

    # Precalculate tf
    tf = math.tan(phif)
    tf2 = tf * tf
    tf4 = tf2 * tf2

    # Precalculate fractional coefficients for x**n in the equations
    # below to simplify the expressions for latitude and longitude.
    x1frac = 1.0 / (Nfpow * cf)

    Nfpow *= Nf  # now equals Nf**2)
    x2frac = tf / (2.0 * Nfpow)

    Nfpow *= Nf  # now equals Nf**3)
    x3frac = 1.0 / (6.0 * Nfpow * cf)

    Nfpow *= Nf  # now equals Nf**4)
    x4frac = tf / (24.0 * Nfpow)

    Nfpow *= Nf  # now equals Nf**5)
    x5frac = 1.0 / (120.0 * Nfpow * cf)

    Nfpow *= Nf  # now equals Nf**6)
    x6frac = tf / (720.0 * Nfpow)

    Nfpow *= Nf  # now equals Nf**7)
    x7frac = 1.0 / (5040.0 * Nfpow * cf)

    Nfpow *= Nf  # now equals Nf**8)
    x8frac = tf / (40320.0 * Nfpow)

    # Precalculate polynomial coefficients for x**n.
    # -- x**1 does not have a polynomial coefficient.
    x2poly = -1.0 - nuf2

    x3poly = -1.0 - 2 * tf2 - nuf2

    x4poly = 5.0 + 3.0 * tf2 + 6.0 * nuf2 - 6.0 * tf2 * nuf2 - 3.0 * (nuf2 * nuf2) - 9.0 * tf2 * (nuf2 * nuf2)

    x5poly = 5.0 + 28.0 * tf2 + 24.0 * tf4 + 6.0 * nuf2 + 8.0 * tf2 * nuf2

    x6poly = -61.0 - 90.0 * tf2 - 45.0 * tf4 - 107.0 * nuf2 + 162.0 * tf2 * nuf2

    x7poly = -61.0 - 662.0 * tf2 - 1320.0 * tf4 - 720.0 * (tf4 * tf2)

    x8poly = 1385.0 + 3633.0 * tf2 + 4095.0 * tf4 + 1575 * (tf4 * tf2)

    # Calculate latitude
    latitude = phif + x2frac * x2poly * (x * x) + x4frac * x4poly * math.pow(x, 4.0) + x6frac * x6poly * math.pow(x, 6.0) + x8frac * x8poly * math.pow(x, 8.0)

    # Calculate longitude
    longitude = lambda0 + x1frac * x + x3frac * x3poly * math.pow(x, 3.0) + x5frac * x5poly * math.pow(x, 5.0) + x7frac * x7poly * math.pow(x, 7.0)

    return latitude, longitude


def utm_xy_to_lat_lon(x, y, zone, northernHemisphere=True):

    x -= 500000.0
    x /= UTMScaleFactor

    # If in southern hemisphere, adjust y accordingly.
    if not northernHemisphere:
        y -= 10000000.0

    y /= UTMScaleFactor

    cmeridian = utm_central_meridian(zone)
    return map_xy_to_lat_lon(x, y, cmeridian)


def transform_entity(entity):
    """
    Parse the entity for properties matching the config and convert the utmToLatLng
    coordinates to LatLong (WSG84)
    """

    global easting_property
    global northing_property
    global zone_property
    global zone_default
    global hemi_property
    global hemi_default
    global hemi_northern_value
    global lat_property
    global long_property
    global include_latlong
    global latlong_property

    #logger.info("easting_property=%s", easting_property)
    #logger.info("northing_property=%s", northing_property)
    #logger.info("zone_property=%s", zone_property)
    #logger.info("zone_default=%s", zone_default)
    #logger.info("hemi_property=%s", hemi_property)
    #logger.info("hemi_default=%s", hemi_default)
    #logger.info("hemi_northern_value=%s", hemi_northern_value)
    #logger.info("lat_property=%s", lat_property)
    #logger.info("long_property=%s", long_property)
    #logger.info("include_latlong=%s", include_latlong)
    #logger.info("latlong_property=%s", latlong_property)

    if easting_property not in entity:
        if logger:
            logger.warning("No easting coordinate found in entity, skipping...")
        return entity

    easting_value = entity.get(easting_property)
    if isinstance(easting_value, list):
        if len(easting_value) > 1:
            if logger:
                logger.warning("Multiple easting coordinates found in entity, skipping...")
            return entity
        easting_value = easting_value[0]
    elif not easting_value:
        if logger:
            logger.warning("skipping due to null easting value...")
        return entity

    if northing_property not in entity:
        if logger:
            logger.warning("No northing coordinate found in entity, skipping...")
        return entity

    northing_value = entity.get(northing_property)
    if isinstance(northing_value, list):
        if len(northing_value) > 1:
            if logger:
                logger.warning("Multiple northing coordinates found in entity, skipping...")
            return entity
        northing_value = northing_value[0]
    elif not northing_value:
        if logger:
            logger.warning("skipping due to null northing value...")
        return entity

    zone_value = entity.get(zone_property, zone_default)
    if isinstance(zone_value, list):
        if len(zone_value) > 1:
            if logger:
                logger.warning("Multiple zone values found in entity, skipping...")
            return entity
        zone_value = zone_value[0]

    hemi_value = entity.get(hemi_property, hemi_default)
    if isinstance(hemi_value, list):
        if len(hemi_value) > 1:
            if logger:
                logger.warning("Multiple hemisphere values found in entity, skipping...")
            return entity
        hemi_value = hemi_value[0]

    # Ready to convert to latlong

    if isinstance(easting_value, str):
        easting_value = easting_value.strip()

    if isinstance(northing_value, str):
        northing_value = northing_value.strip()

    if isinstance(zone_value, str):
        zone_value = zone_value.strip()

    if isinstance(hemi_value, str):
        hemi_value = hemi_value.strip()

    try:
        easting_value = float(easting_value)
    except:
        msg = "Could not convert easting value '%s' to float - format error!" % easting_value
        if logger:
            logger.error(msg)
        raise AssertionError(msg)

    try:
        northing_value = float(northing_value)
    except:
        msg = "Could not convert northing value '%s' to float - format error!" % northing_value
        if logger:
            logger.error(msg)
        raise AssertionError(msg)

    try:
        zone_value = int(zone_value)
    except:
        msg = "Could not convert zone value '%s' to integer - format error!" % zone_value
        if logger:
            logger.error(msg)
        raise AssertionError(msg)

    hemi_value = (hemi_value == hemi_northern_value)

    if logger:
        logger.debug("Converting %s %s, %s %s..." % (easting_value, northing_value, zone_value, hemi_value))

    latitude, longitude = utm_xy_to_lat_lon(easting_value, northing_value, zone_value, hemi_value)

    lat = rad_to_deg(latitude)
    lon = rad_to_deg(longitude)

    if logger:
        logger.debug("Result: %s %s" % (lat, lon))

    entity[lat_property] = lat
    entity[long_property] = lon

    b_include_latlong = (include_latlong.strip().lower() == "true")

    if b_include_latlong:
        entity[latlong_property] = "%s, %s" % (lat, lon)

    return entity


@app.route('/transform', methods=['POST'])
def receiver():
    """ HTTP transform POST handler """

    def generate(entities):
        yield "["
        for index, entity in enumerate(entities):
            if index > 0:
                yield ","

            # Transit decode
            
                
            entity = transform_entity(entity)
            yield entities_to_json(entity)
        yield "]"

    # get entities from request
    req_entities = parse_json_stream(request.stream)

    # Generate the response
    try:
        return Response(generate(req_entities), mimetype='application/json')
    except BaseException as e:
        return Response(status=500, response="An error occured during transform of input")


if __name__ == '__main__':
    
    # Set up logging
    format_string = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logger = logging.getLogger('utmtolatlong-microservice')

    # Log to stdout
    stdout_handler = logging.StreamHandler()
    stdout_handler.setFormatter(logging.Formatter(format_string))
    logger.addHandler(stdout_handler)

    logger.setLevel(logging.DEBUG)

    app.run(debug=True, host='0.0.0.0', port=5001)

