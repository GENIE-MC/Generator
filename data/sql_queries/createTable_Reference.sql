CREATE TABLE IF NOT EXISTS REFERENCE (

  name             VARCHAR(30)   NOT NULL,
  measurement_tag  VARCHAR(30)   DEFAULT '0' NOT NULL,
  iref             INT(2)        DEFAULT  0  NOT NULL,

  authors          VARCHAR(50)  NOT NULL,
  journal          VARCHAR(30)  NOT NULL,
  year             YEAR         NOT NULL,

  PRIMARY KEY (name, measurement_tag, iref)
);
