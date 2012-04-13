CREATE TABLE IF NOT EXISTS MEASUREMENT_HEADER (

  name             VARCHAR(30)   NOT NULL,
  measurement_tag  VARCHAR(30)   DEFAULT '0' NOT NULL,

  observable       VARCHAR(50)   NOT NULL,
  target           VARCHAR(30)   NOT NULL,
  reaction         VARCHAR(30)   NOT NULL,
  A                VARCHAR(5)    DEFAULT '1' NOT NULL,
  exposure         VARCHAR(30),
  exposure_units   VARCHAR(30),
  data_source      VARCHAR(30)   NOT NULL,
  npoints          INT(4)        NOT NULL,
  err_status       INT(1)        DEFAULT '1' NOT NULL,
  nrefs            INT(2)        NOT NULL,

  comment          TEXT,

  PRIMARY KEY (name, measurement_tag)
);
