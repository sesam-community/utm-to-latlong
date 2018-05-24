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


def utm_to_lat_long(zone, easting, northing, northernHemisphere=True):
    if not northernHemisphere:
        northing = 10000000 - northing

    a = 6378137
    e = 0.081819191
    e1sq = 0.006739497
    k0 = 0.9996

    arc = northing / k0
    mu = arc / (a * (1 - math.pow(e, 2) / 4.0 - 3 * math.pow(e, 4) / 64.0 - 5 * math.pow(e, 6) / 256.0))

    ei = (1 - math.pow((1 - e * e), (1 / 2.0))) / (1 + math.pow((1 - e * e), (1 / 2.0)))

    ca = 3 * ei / 2 - 27 * math.pow(ei, 3) / 32.0

    cb = 21 * math.pow(ei, 2) / 16 - 55 * math.pow(ei, 4) / 32
    cc = 151 * math.pow(ei, 3) / 96
    cd = 1097 * math.pow(ei, 4) / 512
    phi1 = mu + ca * math.sin(2 * mu) + cb * math.sin(4 * mu) + cc * math.sin(6 * mu) + cd * math.sin(8 * mu)

    n0 = a / math.pow((1 - math.pow((e * math.sin(phi1)), 2)), (1 / 2.0))

    r0 = a * (1 - e * e) / math.pow((1 - math.pow((e * math.sin(phi1)), 2)), (3 / 2.0))
    fact1 = n0 * math.tan(phi1) / r0

    _a1 = 500000 - easting
    dd0 = _a1 / (n0 * k0)
    fact2 = dd0 * dd0 / 2

    t0 = math.pow(math.tan(phi1), 2)
    Q0 = e1sq * math.pow(math.cos(phi1), 2)
    fact3 = (5 + 3 * t0 + 10 * Q0 - 4 * Q0 * Q0 - 9 * e1sq) * math.pow(dd0, 4) / 24

    fact4 = (61 + 90 * t0 + 298 * Q0 + 45 * t0 * t0 - 252 * e1sq - 3 * Q0 * Q0) * math.pow(dd0, 6) / 720

    lof1 = _a1 / (n0 * k0)
    lof2 = (1 + 2 * t0 + Q0) * math.pow(dd0, 3) / 6.0
    lof3 = (5 - 2 * Q0 + 28 * t0 - 3 * math.pow(Q0, 2) + 8 * e1sq + 24 * math.pow(t0, 2)) * math.pow(dd0, 5) / 120
    _a2 = (lof1 - lof2 + lof3) / math.cos(phi1)
    _a3 = _a2 * 180 / math.pi

    latitude = 180 * (phi1 - fact1 * (fact2 + fact3 + fact4)) / math.pi

    if not northernHemisphere:
        latitude = -latitude

    longitude = ((zone > 0) and (6 * zone - 183.0) or 3.0) - _a3

    return (latitude, longitude)



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

    lat, lon = utm_to_lat_long(zone_value, easting_value, northing_value, hemi_value)

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

