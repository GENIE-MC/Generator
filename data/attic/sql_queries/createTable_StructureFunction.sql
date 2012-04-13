CREATE TABLE IF NOT EXISTS STRUCTURE_FUNCTION (

  name             VARCHAR(30)    NOT NULL,
  measurement_tag  VARCHAR(30)    DEFAULT '0' NOT NULL,

  sf               FLOAT          NOT NULL,
  p                VARCHAR(20)    NOT NULL,
  R                VARCHAR(10)    DEFAULT 'QCD' NOT NULL,
  x                FLOAT          NOT NULL,
  Q2               FLOAT          NOT NULL,
  stat_err_p       FLOAT,
  stat_err_m       FLOAT,
  syst_err_p       FLOAT,
  syst_err_m       FLOAT,

  PRIMARY KEY (name, measurement_tag, x, Q2, p, R)
);
