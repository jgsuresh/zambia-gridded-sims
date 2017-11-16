from datetime import *
from dateutil import relativedelta

def convert_to_day(convert_date, ref_date, date_format = "%Y-%m-%d"):
    # Converts date to day of simulation starting from reference date
    # Uses actual calendar dates
    date_delta = datetime.strptime(convert_date, date_format) - \
                 datetime.strptime(ref_date, date_format)

    return date_delta.days

def convert_to_date(convert_day, ref_date, date_format = "%Y-%m-%d"):
    # Converts day of simulation starting from reference date into date
    # Uses actual calendar dates
    full_date = datetime.strptime(ref_date, date_format) + relativedelta.relativedelta(days=int(convert_day))

    return datetime.strftime(full_date, date_format)

def convert_to_day_365(convert_date, ref_date, date_format = "%Y-%m-%d"):
    # Converts date to day of simulation starting from reference date
    # Assumes a calendar year has exactly 365 days
    return 365*(datetime.strptime(convert_date, date_format).year - \
                datetime.strptime(ref_date, date_format).year) + \
           datetime.strptime(convert_date, date_format).timetuple().tm_yday

def convert_to_date_365(convert_day, ref_date, date_format = "%Y-%m-%d"):
    # Converts day of simulation starting from reference date into date
    # Assumes a calendar year has exactly 365 days

    convert_year = (int(convert_day)/365) + datetime.strptime(ref_date, date_format).year
    convert_day = convert_day-(int(convert_day)/365)*365
    return datetime.strftime(datetime(convert_year, 1, 1) + timedelta(convert_day - 1), date_format)
