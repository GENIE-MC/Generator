CREATE TABLE IF NOT EXISTS EXP_INFO (

  name             VARCHAR(30)   NOT NULL,

  comment          TEXT,

  facility         VARCHAR(30)   NOT NULL,
  detector         VARCHAR(30)   NOT NULL,
  beam             VARCHAR(30)   NOT NULL,
  target           VARCHAR(30)   NOT NULL,

  year_start       YEAR,
  year_end         YEAR,

  exposure         VARCHAR(30),
  exposure_units   VARCHAR(30),

  energy_min       FLOAT,
  energy_max       FLOAT,
  energy_units     VARCHAR(30)   DEFAULT 'GeV',
  energy_frame     VARCHAR(30)   DEFAULT 'lab',

  PRIMARY KEY (name)
);
